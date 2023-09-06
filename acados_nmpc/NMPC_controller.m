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
        W_x = diag([100, 1, 1,  1e-3]);
        W_x_e = 2*diag([100, 1, 1, 1e-3]);
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
        %         ocp_model;
        %         ocp_opts;
        ocp_solver;
        cost_function_vect;

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

        % Delay compensation parameters
        delay_compensation; % Delay to compensate with the controller [s]
        delay_buff_comp;
        u_buff_contr;

        plant;

        v_alpha;
        d_v_bound;
        t_angle0;

    end

    methods

        % Constructor
        function self = NMPC_controller(name,plant, sample_time, Hp)
            self@casadi.Callback();

            % Acados initialitation
            %self.init();
            check_acados_requirements()

            % Save the plant model
            self.sym_model = plant.sym_model;
            self.sym_model.name = plant.name;
            self.initial_condition = zeros(plant.sym_model.nx,1);

            % Constraints
%             self.h_constr_ub = [0.9*plant.slider_params.ywidth/2 self.u_n_ub self.u_t_ub];
%             self.h_constr_lb = [-0.9*plant.slider_params.ywidth/2 self.u_n_lb self.u_t_lb];
            self.h_constr_ub = [10 self.u_n_ub self.u_t_ub];
            self.h_constr_lb = [-10 self.u_n_lb self.u_t_lb];

            % Controller parameters
            self.Hp = Hp;
            self.sample_time = sample_time;
            self.T = self.Hp*self.sample_time;

            self.plant = plant;

%             self.set_v_alpha(0.035);%0.005*200);
%             self.d_v_bound = 0.0045;

%             self.set_v_alpha(0.005*200);
%             self.d_v_bound = 0;
            self.set_v_alpha(1*0.5719);%0.005*200);
            self.d_v_bound = 0.0;
            self.t_angle0 = 3+0*2.831;

            construct(self, name);
        end

        function set_delay_comp(self,delay)
            self.delay_compensation = delay;
            self.delay_buff_comp = ceil(self.delay_compensation/self.sample_time);
            self.u_buff_contr = zeros*ones(self.sym_model.nu, self.delay_buff_comp);
        end

        function xk_sim = delay_buffer_sim(self, plant, x)
            xk_sim = x;
            for k = 1 : self.delay_buff_comp
%                 x_dot_sim = plant.eval_model(xk_sim,self.u_buff_contr(:,end-k+1));
                x_dot_sim = plant.evalModelVariableShape(xk_sim,self.u_buff_contr(:,end-k+1));
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

%             self.ocp_solver.set('constr_lbu', [self.u_n_lb; self.u_t_lb]); % lower bound on h
%             self.ocp_solver.set('constr_ubu', [self.u_n_ub; self.u_t_ub]);  % upper bound on h

%             self.ocp_solver.set('constr_lbu', self.u_n_lb); % lower bound on h
%             self.ocp_solver.set('constr_ubu', self.u_n_ub);  % upper bound on h


            self.ocp_solver.set('constr_lh', [self.h_constr_lb(2:end-1), 2*self.u_t_lb, -0]); % lower bound on h
            self.ocp_solver.set('constr_uh', [self.h_constr_ub(2:end-1), 0, 2*self.u_t_ub]);  % upper bound on h

%             self.ocp_solver.set('constr_lh', [self.h_constr_lb(2:end)]);% 2*self.u_t_lb 0]); % lower bound on h
%             self.ocp_solver.set('constr_uh', [self.h_constr_ub(2:end)]);% 0 2*self.u_t_ub]);  % upper bound on h
        end
        
        
        function clear_variables(self)
            % Usefull method to clean some variables for the new experiment
            self.utraj = [];
            self.xtraj = [];
            self.ptraj = [];
            self.y_ref = [];
            self.cost_function_vect = [];
        end

        function update_cost_function(self,W_x,W_u,W_x_e,initial_step, final_step)
            self.ocp_solver.set('cost_W', W_x_e,self.Hp);

            for i = initial_step : final_step
                self.ocp_solver.set('cost_W', blkdiag(W_x,W_u),i);
            end
            self.W_x = W_x;
            self.W_u = W_u;
            self.W_x_e=W_x_e;


        end

        function initial_condition_update(self,new_initial_condition)
            % Update the initial condition for a new control
            % Input: x0
            self.initial_condition = new_initial_condition;
            self.ocp_solver.set('constr_x0', new_initial_condition); %%%%%% MC: perchè si setta su ocp_model e non su ocp_solver?
            self.clear_variables();
        end

        function ocp_model = create_ocp_model(self)
            import casadi.*
            % acados ocp model
            ocp_model = acados_ocp_model();

            ocp_model.set('name', self.sym_model.name);

            % symbolics
            ocp_model.set('sym_x', self.sym_model.sym_x);
            ocp_model.set('sym_u', self.sym_model.sym_u);

            ocp_model.set('cost_type', 'linear_ls');
            ocp_model.set('cost_type_e', 'linear_ls');

            Vx = zeros(self.sym_model.nx+self.sym_model.nu, self.sym_model.nx);
            Vx_e = zeros(self.sym_model.nx, self.sym_model.nx);
            Vu = zeros(self.sym_model.nx+self.sym_model.nu, self.sym_model.nu);
            for ii=1:self.sym_model.nx
                Vx(ii,ii)=1.0;
            end
            for ii=1:self.sym_model.nu
                Vu(self.sym_model.nx+ii,ii)=1.0;
            end
            for ii=1:self.sym_model.nx
                Vx_e(ii,ii)=1.0;
            end

            %             ocp_model.set('cost_Vx',[eye(self.sym_model.nx); zeros(self.sym_model.nu,self.sym_model.nx)]);
            %             ocp_model.set('cost_Vu',[zeros(self.sym_model.nx,self.sym_model.nu); eye(self.sym_model.nu)]);
            ocp_model.set('cost_Vx',Vx);
            ocp_model.set('cost_Vu',Vu);
            ocp_model.set('cost_Vz',zeros(self.sym_model.nx+self.sym_model.nu,0));

            % var
            ocp_model.set('cost_W',blkdiag(self.W_x,self.W_u));
            ocp_model.set('cost_y_ref',zeros(self.sym_model.nx+self.sym_model.nu,1));
            %             ocp_model.set('cost_y_ref',zeros(self.sym_model.nx,1));

            %             ocp_model.set('cost_Vx_e',eye(self.sym_model.nx));
            ocp_model.set('cost_Vx_e',Vx_e);

            % var
            %%%%%% CAMBIO PER W_E
            ocp_model.set('cost_W_e',(self.W_x_e));
            ocp_model.set('cost_y_ref_e',zeros(self.sym_model.nx,1));


            ocp_model.set('T', self.T);


            ocp_model.set('dyn_type', 'explicit');
            ocp_model.set('dyn_expr_f', self.sym_model.expr_f_expl);

            % constraints
%             expr_h = vertcat(self.sym_model.sym_x(end),self.sym_model.sym_u);  %varables to be constrained
            s_mod = Function('s_mod',{self.sym_model.sym_x(4)},{(self.sym_model.sym_x(4)<0)*(self.plant.SP.b)+mod(self.sym_model.sym_x(4),self.plant.SP.b)});
            t_angle_fun = self.plant.SP.FC_angle_dot(1);
            t_angle = Function('t_angle',{self.sym_model.sym_x(4)},{abs(t_angle_fun(s_mod(self.sym_model.sym_x(4))))});
            v_bound = min((self.v_alpha/(abs(t_angle(s_mod(self.sym_model.sym_x(4)))-self.t_angle0) + 0.0001)) + self.d_v_bound,self.u_t_ub);
            v_bound_f = Function('v_bound_f',{self.sym_model.sym_x(4)},{v_bound});
            u_t_bound_n = Function('u_t_bound_n',{self.sym_model.sym_x(4)},{(self.sym_model.sym_u(2)-v_bound_f(s_mod(self.sym_model.sym_x(4))))});
            u_t_bound_p = Function('u_t_bound_p',{self.sym_model.sym_x(4)},{(self.sym_model.sym_u(2)+v_bound_f(s_mod(self.sym_model.sym_x(4))))});
% 
%             expr_h = vertcat(self.sym_model.sym_u(1),self.sym_model.sym_u(2)); % (self.sym_model.sym_u(2)+v_bound_f(s_mod(self.sym_model.sym_x(4)))));
            expr_h = vertcat(self.sym_model.sym_u(1),u_t_bound_n(self.sym_model.sym_x(4)),u_t_bound_p(self.sym_model.sym_x(4))); % (self.sym_model.sym_u(2)+v_bound_f(s_mod(self.sym_model.sym_x(4)))));
    
            ocp_model.set('constr_type', 'bgh');
%             ocp_model.set('constr_type', 'auto');
            ocp_model.set('constr_expr_h', expr_h);
% 
%             ocp_model.set('constr_lh', self.h_constr_lb); % lower bound on h
%             ocp_model.set('constr_uh', self.h_constr_ub);  % upper bound on h
% 
            ocp_model.set('constr_lh', [self.h_constr_lb(2:end-1) 2*self.u_t_lb -0]); % lower bound on h
            ocp_model.set('constr_uh', [self.h_constr_ub(2:end-1) 0 2*self.u_t_ub]);  % upper bound on h


%             ocp_model.set('constr_lh', [self.h_constr_lb(2:end)]);% 2*self.u_t_lb 0]); % lower bound on h
%             ocp_model.set('constr_uh', [self.h_constr_ub(2:end)]);% 0 2*self.u_t_ub]);  % upper bound on h


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
            %             self.ocp_model = ocp_model;
            % ... see ocp_model.model_struct to see what other fields can be set
        end

        function ocp_opts = create_ocp_opts(self)
            field_s = ["nlp_solver", "qp_solver", "sim_method",  "globalization", "codgen_model", "compile_model", "compile_interface"];
            values_s = ["sqp", "partial_condensing_hpipm", "erk", "merit_backtracking", "true", "false", "false"];
            %solver_params_s = dictionary(field_s,values_s);

            field_d = ["qp_solver_cond_N", "nlp_solver_max_iter","line_search_use_sufficient_descent","nlp_solver_tol_stat","nlp_solver_tol_eq","nlp_solver_tol_ineq","nlp_solver_tol_comp"];
            values_d = [5, 30, 1,1e-6,1e-6,1e-6,1e-6];
            %solver_params_d = dictionary(field_d,values_d);

            % acados ocp set opts
            ocp_opts = acados_ocp_opts();
            ocp_opts.set('param_scheme_N', self.Hp);
            %ent = solver_params_s.keys;
            %val = solver_params_s.values;
            for index = 1:length(field_s)
                ocp_opts.set(field_s(index),values_s(index));
            end
 
%             ent = solver_params_d.keys;
%             val = solver_params_d.values;
            for index = 1:length(field_d)
                ocp_opts.set(field_d(index),values_d(index));
            end


            %ocp_opts.set('codgen_model','false');
            ocp_opts.set('compile_interface','auto');
            ocp_opts.set('output_dir',fullfile(pwd,'build'));
            %             self.ocp_opts = ocp_opts;
            % ... see ocp_opts.opts_struct to see what other fields can be set
        end

        function create_ocp_solver(self)
            %create the ocp solver
            self.ocp_solver = acados_ocp(self.create_ocp_model(), self.create_ocp_opts());
        end

        function y_ref_k = get_y_ref(self, index_ref)
            if(index_ref>length(self.y_ref))
                y_ref_k = self.y_ref(:,end);
            else
                y_ref_k = self.y_ref(:,index_ref);
            end
        end

        function set_v_alpha(self,alpha)
            self.v_alpha = alpha;
        end

        function [v_bound, t_angle] = update_tangential_velocity_bounds(self,s)
            s = mod(s,self.plant.SP.b);
            t_angle = abs(self.plant.SP.getAngleCurvatures(s));
            v_bound = min((self.v_alpha/(abs(t_angle-self.t_angle0)+0.0001)) + self.d_v_bound,self.u_t_ub);
            

%             self.ocp_solver.set('constr_lbu', [self.u_n_lb; -v_bound]); % lower bound on h
%             self.ocp_solver.set('constr_ubu', [self.u_n_ub; v_bound]);  % upper bound on h
        end
        
        function u = solve(self,x0, index_time)
            % update initial state
            %             tic
            x0(4) = mod(x0(4),self.plant.SP.b);

            self.ocp_solver.set('constr_x0', x0);
            
%             alpha = self.plant.SP.getNormalizedCurvature(x0(4));
%             scale_alpha = 0.7;
%             self.ocp_solver.set('constr_lbu', [self.u_n_lb; (1-scale_alpha*alpha)*self.u_t_lb]); % lower bound on h
%             self.ocp_solver.set('constr_ubu', [self.u_n_ub; (1-scale_alpha*alpha)*self.u_t_ub]);  % upper bound on h
%             self.update_tangential_velocity_bounds(x0(4));

            % reference
            for k=0:self.Hp-1
                y_ref_k = self.get_y_ref(index_time+k);
                self.ocp_solver.set('cost_y_ref', y_ref_k, k);  %reference for middle samples of the prediction horizon
            end

            self.ocp_solver.set('cost_y_ref_e', y_ref_k(1:self.sym_model.nx), self.Hp);

            % set initial guess
            if(isempty(self.utraj) || isempty(self.xtraj))
                self.xtraj = zeros(self.sym_model.nx, self.Hp+1);
                self.utraj = repmat([self.u_n_lb;0],1,self.Hp);
                self.ptraj = zeros(self.sym_model.nx, self.Hp);
            end

            v_bounds = self.update_tangential_velocity_bounds(x0(4));
            if abs(self.utraj(2,1)) > v_bounds 
                disp("Violated constraint")
                
                ut_old = self.utraj(2,1);
                self.utraj(2,1) = sign(ut_old)*v_bounds;
                self.utraj(1,1) = self.utraj(2,1)*self.utraj(1,1)/ut_old;
            end

            self.xtraj(:,1) = x0;
            for i_x = 2 : size(self.xtraj,2)
%                 self.xtraj(:,i_x) = self.plant.eval_model(self.xtraj(:,i_x-1),self.utraj(:,i_x-1));
                xtraj_dot = self.plant.evalModelVariableShape(self.xtraj(:,i_x-1),self.utraj(:,i_x-1));
                self.xtraj(:,i_x) = self.xtraj(:,i_x-1) + self.sample_time*xtraj_dot; 
                v_bounds = self.update_tangential_velocity_bounds(self.xtraj(4,i_x));
                if i_x == size(self.xtraj,2)
                    break;
                end
                if abs(self.utraj(2,i_x)) > v_bounds
                    ut_old = self.utraj(2,i_x);
                    self.utraj(2,i_x) = sign(ut_old)*v_bounds;
                    self.utraj(1,i_x) = self.utraj(2,i_x)*self.utraj(1,i_x)/ut_old;
                end
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


            self.utraj = [self.utraj(:,2:end) self.utraj(:,end)];
            self.xtraj = [self.xtraj(:,2:end) self.xtraj(:,end)];
            self.ptraj = [self.ptraj(:,2:end) self.ptraj(:,end)];

            % status = ocp.get('status'); % 0 - success
            % ocp.print('stat')
            u = self.ocp_solver.get('u', 0);
            [v_bounds, t_angle] = self.update_tangential_velocity_bounds(x0(4));

%             if abs(u(2)) > v_bounds
%                 ut_old = u(2);
%                 u(2) = sign(ut_old)*v_bounds;
%                 u(1) = (u(2)/ut_old)*u(1);
            if t_angle > 120
                u(1) = -0.001;
            end
%                 
%             end


            self.cost_function_vect = [self.cost_function_vect; self.ocp_solver.get_cost];

            %             toc
        end

        function set_reference_trajectory(self,y_ref)
            % Update reference trajectory
            self.y_ref = [repmat([self.initial_condition; 0; 0],1,self.delay_buff_comp) y_ref];%[zeros(size(y_ref,1),self.delay_buff_comp) y_ref];
            self.y_ref(self.sym_model.nx+self.sym_model.nu,1:self.delay_buff_comp) = self.y_ref(self.sym_model.nx+self.sym_model.nu,self.delay_buff_comp+1);
            %             self.y_ref = y_ref(:,self.delay_buff_comp:end);
        end


    end
end