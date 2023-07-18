classdef NMPC_controller < casadi.Callback
    % This class creates a nonlinear mpc controller by using acados.
    %
    % NMPC_controller Properties:
    %    sample_time - Sample time of the controller
    %    W_x - state weight matrix
    %    W_u - control weight matrix
    %
    % NMPC_controller Methods:
    %    NMPC_controller - Constructor
    %

    properties

        % Weight matrices for objective function
        W_x = diag([100, 1, 1, 1e-3, 1e-3]);
        W_u = diag([1 1]);

        % Costraints
        h_constr_lb = [];  % lower bound constraint on the constrained variables h
        h_constr_ub = [];  % upper bound constraint on the constrained variables h
        u_n_ub = 0.05;
        u_t_ub = 0.05;
        u_n_lb = 0;
        u_t_lb = -0.05;

        % Solver parameters
        solver_params = struct;
        ocp_model;
        ocp_opts;
        ocp_solver;

        % Symbolic plant model
        sym_model = struct;
        initial_condition;

        % Controller params
        Hp;                     % prediction horizon
        sample_time;            % controller sample time [s]
        T;                      % time horizon length [s]

        % Previous solutions
        utraj;
        xtraj;
        ptraj;

        % Reference
        y_ref;             % Trajectory points
        index_ref = 0;     % utility variable for the trajectory tracking

        % Delay compensation parameters
        delay_compensation; % Delay to compensate with the controller [s]
        delay_buff_comp;
        u_buff_contr;

    end

    methods(Static)

        function solver_params = set_solver_params(self)
            % Solver parameters
            solver_params.nlp_solver = 'sqp_rti';                      % sqp, sqp_rti
            solver_params.qp_solver = 'partial_condensing_hpipm';  % full_condensing_hpipm, partial_condensing_hpipm, full_condensing_qpoases, full_condensing_daqp
            solver_params.qp_solver_cond_N = 5;                    % for partial condensing
            solver_params.sim_method = 'erk';                      % integrator type : erk, irk, irk_gnsf
            self.solver_params = solver_params;
        end
    end
    methods

        % Constructor
        function self = NMPC_controller(name,plant, sample_time, Hp)
            self@casadi.Callback();

            % Acados initialitation
            %self.init();
            check_acados_requirements()

            % Set solver parameters
            self.solver_params = self.set_solver_params();

            % Save the plant model
            self.sym_model = plant.sym_model;
            self.sym_model.name = plant.name;
            self.initial_condition = zeros(plant.sym_model.nx,1);

            % Constraints
            self.h_constr_ub = [0.9*plant.slider_params.ywidth/2 self.u_n_ub self.u_t_ub];
            self.h_constr_lb = [-0.9*plant.slider_params.ywidth/2 self.u_n_lb self.u_t_lb];

            % Controller parameters
            self.Hp = Hp;
            self.sample_time = sample_time;
            self.T = self.Hp*self.sample_time;

            construct(self, name);
        end

        function set_delay_comp(self,delay)
            self.delay_compensation = delay;
            self.delay_buff_comp = ceil(self.delay_compensation/self.sample_time);
            self.u_buff_contr = zeros(self.sym_model.nu, self.delay_buff_comp);
        end

        function xk_sim = delay_buffer_sim(self, plant, x)
            xk_sim = x;
            for k = 1 : self.delay_buff_comp
                x_dot_sim = plant.eval_model(xk_sim,self.u_buff_contr(:,end-k+1));
                x_sim = xk_sim + self.sample_time*x_dot_sim;
                xk_sim = x_sim;
            end
        end

        function update_constraints(self, u_n_ub, u_t_ub, u_n_lb, u_t_lb)
            self.u_n_lb = u_n_lb;
            self.u_n_ub = u_n_ub;
            self.u_t_lb = u_t_lb;
            self.u_t_ub = u_t_ub;
            self.h_constr_ub = [self.h_constr_ub(1) self.u_n_ub self.u_t_ub];
            self.h_constr_lb = [self.h_constr_lb(1) self.u_n_lb self.u_t_lb];

            self.ocp_model.set('constr_lh', self.h_constr_lb); % lower bound on h
            self.ocp_model.set('constr_uh', self.h_constr_ub);  % upper bound on h
            self.ocp_solver.set('constr_lbu', [self.u_n_lb; self.u_t_lb]); % lower bound on h
            self.ocp_solver.set('constr_ubu', [self.u_n_ub; self.u_t_ub]);  % upper bound on h
        end

        function clear_variables(self)
            % Usefull method to clean some variables for the new experiment
            self.utraj = [];
            self.xtraj = [];
            self.ptraj = [];
            self.index_ref = 0;
            self.y_ref = [];
        end

        function update_cost_function(self,W_x,W_u,W_x_e,initial_step, final_step)

            if(initial_step == self.Hp)
                % Last step of the prediction horizon, means W_e
                self.ocp_solver.set('cost_W', W_x_e,self.Hp);
            else
                for i = initial_step : final_step
                    self.ocp_solver.set('cost_W', blkdiag(W_x,W_u),i);
                end
            end
        end

        function initial_condition_update(self,new_initial_condition)
            % Update the initial condition for a new control
            % Input: x0
            self.initial_condition = new_initial_condition;
            self.ocp_model.set('constr_x0', new_initial_condition);
            self.clear_variables();
        end

        function ocp_model = create_ocp_model(self)
            % acados ocp model
            ocp_model = acados_ocp_model();

            ocp_model.set('name', self.sym_model.name);
            ocp_model.set('T', self.T);

            % symbolics
            ocp_model.set('sym_x', self.sym_model.sym_x);
            ocp_model.set('sym_u', self.sym_model.sym_u);
            ocp_model.set('sym_xdot', self.sym_model.sym_xdot);

            % cost
            expr_ext_cost_e = self.sym_model.sym_x'* self.W_x * self.sym_model.sym_x;
            expr_ext_cost = expr_ext_cost_e + self.sym_model.sym_u' * self.W_u * self.sym_model.sym_u;

            % nonlinear least sqares
            %cost_expr_y = vertcat(self.sym_model.sym_x, self.sym_model.sym_u);
            %W = blkdiag(self.W_x, self.W_u);
            self.sym_model.cost_expr_y_e = self.sym_model.sym_x;
            self.sym_model.W_e = self.W_x;

            ocp_model.set('cost_expr_ext_cost', expr_ext_cost);
            ocp_model.set('cost_expr_ext_cost_e', expr_ext_cost_e);

            ocp_model.set('dyn_type', 'explicit');
            ocp_model.set('dyn_expr_f', self.sym_model.expr_f_expl);

            % constraints
            expr_h = vertcat(self.sym_model.sym_x(end),self.sym_model.sym_u);  %varables to be constrained

            ocp_model.set('constr_type', 'auto');
            ocp_model.set('constr_expr_h', expr_h);

            ocp_model.set('constr_lh', self.h_constr_lb); % lower bound on h
            ocp_model.set('constr_uh', self.h_constr_ub);  % upper bound on h

            % % on the state
            % ocp_model.set('constr_Jbx',eye(nx));
            % ocp_model.set('constr_lbx',[-100 -100 -100 -100 -model.slider.ywidth/2]);
            % ocp_model.set('constr_ubx',[100 100 100 100 model.slider.ywidth/2]);
            %
            % % inputs
            % ocp_model.set('constr_Jbu',eye(nu));
            % ocp_model.set('constr_lbu',U_min);
            % ocp_model.set('constr_ubu',U_max);

            ocp_model.set('constr_x0', self.initial_condition);
            self.ocp_model = ocp_model;
            % ... see ocp_model.model_struct to see what other fields can be set
        end

        function ocp_opts = create_ocp_opts(self)
            % acados ocp set opts
            ocp_opts = acados_ocp_opts();
            ocp_opts.set('param_scheme_N', self.Hp);
            ocp_opts.set('nlp_solver', self.solver_params.nlp_solver);
            ocp_opts.set('sim_method', self.solver_params.sim_method);
            ocp_opts.set('qp_solver', self.solver_params.qp_solver);
            ocp_opts.set('qp_solver_cond_N', self.solver_params.qp_solver_cond_N);
            ocp_opts.set('ext_fun_compile_flags', ''); % '-O2'
            ocp_opts.set('compile_model','false');
            %ocp_opts.set('codgen_model','false');
            ocp_opts.set('compile_interface','false');
            ocp_opts.set('output_dir',fullfile(pwd,'build'));
            self.ocp_opts = ocp_opts;
            % ... see ocp_opts.opts_struct to see what other fields can be set
        end

        function create_ocp_solver(self)
            %create the ocp solver
            self.ocp_solver = acados_ocp(self.create_ocp_model(), self.create_ocp_opts());
        end
        
        function u = solve(self,x0)
            % update initial state
%             tic
            self.ocp_solver.set('constr_x0', x0);

            if (self.Hp + self.index_ref) ~= length(self.y_ref)
                % reference
                for k=0:self.Hp-1
                    self.ocp_solver.set('cost_y_ref', self.y_ref(:,self.index_ref+k+1), k);  %reference for middle samples of the prediction horizon
                end
                self.index_ref = self.index_ref + 1;
                self.ocp_solver.set('cost_y_ref_e', self.y_ref(1:self.sym_model.nx,self.index_ref+k+1), self.Hp); % desired state at the end of the prediction horizon
            end
            self.ocp_solver.set('cost_y_ref_e', self.y_ref(1:self.sym_model.nx,self.index_ref+self.Hp), self.Hp);

            % set initial guess
            if(isempty(self.utraj) || isempty(self.xtraj))
                self.xtraj = zeros(self.sym_model.nx, self.Hp+1);
                self.utraj = zeros(self.sym_model.nu, self.Hp);
                self.ptraj = ones(self.sym_model.nx, self.Hp);
            end
            self.ocp_solver.set('init_x', self.xtraj);
            self.ocp_solver.set('init_u', self.utraj);
            self.ocp_solver.set('init_pi', self.ptraj);

            %self.ocp_solver.set('init_pi', ones(self.sym_model.nx, self.Hp));
            %self.ocp_solver.set('constr_lbx', x0, 0)

            self.ocp_solver.solve();    %solve the optimum problem

            % get solution
            self.utraj = self.ocp_solver.get('u');
            self.xtraj = self.ocp_solver.get('x');
            self.ptraj = self.ocp_solver.get('pi');
%             self.utraj = [self.utraj(:,2:end) self.utraj(:,1)];
%             self.xtraj = [self.xtraj(:,2:end) self.xtraj(:,1)];
%             self.ptraj = [self.ptraj(:,2:end) self.ptraj(:,1)];

            % status = ocp.get('status'); % 0 - success
            % ocp.print('stat')
            u = self.ocp_solver.get('u', 0);

%             toc
        end
        
        function set_reference_trajectory(self,y_ref)
            % Update reference trajectory
            self.y_ref = [zeros(size(y_ref,1),self.delay_buff_comp) y_ref];
            self.y_ref(4,1:self.delay_buff_comp) = self.y_ref(4,self.delay_buff_comp+1);
%             self.y_ref = y_ref(:,self.delay_buff_comp:end);
        end


    end
end