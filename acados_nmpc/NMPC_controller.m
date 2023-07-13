classdef NMPC_controller < casadi.Callback

    properties
        % Weights
        W_x = diag([100, 1, 1, 1e-3, 1e-3]);
        W_u = diag([1 1]);

        % Costraints
        h_constr_lb = [];  % lower bound constraint on the constrained variables h
        h_constr_ub = [];  % upper bound constraint on the constrained variables h

        % Solver parameters
        solver_params = struct;
        ocp_model;
        ocp_opts;
        ocp_solver;

        % Symbolic plant model
        sym_model = struct;
        initial_condition;

        % Controller params
        N = 20;                 %prediction horizon
        %sample_time = 0.050;    %controller sample time [s]
        T = 1;                  % time horizon length [s]

        % Previous solutions
        utraj;
        xtraj;

        % Reference
        y_ref;             % Trajectory points
        index_ref = 0;     % utility variable for the trajectory tracking

    end

    methods(Static)
        %         function init()
        %         % Acados initialitation
        %             env_vars_acados
        %             check_acados_requirements()
        %         end
        function solver_params = set_solver_params(self)
            % Solver parameters
            solver_params.nlp_solver = 'sqp';                      % sqp, sqp_rti
            solver_params.qp_solver = 'partial_condensing_hpipm';  % full_condensing_hpipm, partial_condensing_hpipm, full_condensing_qpoases, full_condensing_daqp
            solver_params.qp_solver_cond_N = 5;                    % for partial condensing
            solver_params.sim_method = 'erk';                      % integrator type : erk, irk, irk_gnsf
            self.solver_params = solver_params;
        end
    end
    methods

        % Constructor
        function self = NMPC_controller(name,plant,linux_set)
            self@casadi.Callback();
            % Acados initialitation
            %self.init();
            if linux_set == 0
                env_vars_acados
            end
            check_acados_requirements()

            % Set solver parameters
            self.solver_params = self.set_solver_params();

            % Save the plant model
            self.sym_model = plant.sym_model;
            self.sym_model.name = plant.name;
            self.initial_condition = zeros(plant.sym_model.nx,1);

            % Constraints
            self.h_constr_ub = [plant.slider_params.ywidth/2 0.05 0.05];
            self.h_constr_lb = [-plant.slider_params.ywidth/2 0 -0.05];

            construct(self, name);
        end

                function clear_variables(self)
            % clean some variables for the new experiment
            self.utraj = [];
            self.xtraj = [];
            self.index_ref = 0;
            self.y_ref = [];
        end

        function update_cost_function(self,W_x,W_u)
            % Cost function update
            self.W_x = W_x;
            self.W_u = W_u;

            expr_ext_cost_e = self.sym_model.sym_x'* self.W_x * self.sym_model.sym_x;
            expr_ext_cost = expr_ext_cost_e + self.sym_model.sym_u' * self.W_u * self.sym_model.sym_u;

            self.ocp_model.set('cost_expr_ext_cost', expr_ext_cost);
            self.ocp_model.set('cost_expr_ext_cost_e', expr_ext_cost_e);
        end
        function initial_condition_update(self,new_initial_condition)
            % Update the initial condition for a new control
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
            ocp_opts.set('param_scheme_N', self.N);
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
            tic
            self.ocp_solver.set('constr_x0', x0);
           
            if (self.N + self.index_ref) ~= length(self.y_ref)                    
                % reference            
                for k=0:self.N-1  
                   self.ocp_solver.set('cost_y_ref', self.y_ref(:,self.index_ref+k+1), k);  %reference for middle samples of the prediction horizon 
                end
                self.index_ref = self.index_ref + 1;
                self.ocp_solver.set('cost_y_ref_e', self.y_ref(1:self.sym_model.nx,self.index_ref+k+1), self.N); % desired state at the end of the prediction horizon
            end 
            self.ocp_solver.set('cost_y_ref_e', self.y_ref(1:self.sym_model.nx,self.index_ref+self.N), self.N);
            
            % set initial guess
            if(isempty(self.utraj) || isempty(self.xtraj))
                self.xtraj = zeros(self.sym_model.nx, self.N+1);
                self.utraj = zeros(self.sym_model.nu, self.N);
            end 
            self.ocp_solver.set('init_x', self.xtraj);
            self.ocp_solver.set('init_u', self.utraj);
            
            %self.ocp_solver.set('init_pi', ones(self.sym_model.nx, self.N)); 
            %self.ocp_solver.set('constr_lbx', x0, 0) 

            self.ocp_solver.solve();    %solve the optimum problem

            % get solution
            self.utraj = self.ocp_solver.get('u');
            self.xtraj = self.ocp_solver.get('x');

            % status = ocp.get('status'); % 0 - success
            % ocp.print('stat')
            u = self.ocp_solver.get('u', 0);
            toc
        end
        function set_reference_trajectory(self,y_ref)
            % Update reference trajectory
            self.y_ref = y_ref;
        end 


    end
end