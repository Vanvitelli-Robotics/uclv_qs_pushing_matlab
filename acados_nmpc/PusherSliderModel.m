classdef PusherSliderModel < casadi.Callback
    
    properties
      
        % Pusher_Slider parameters
        slider_params = struct;

        % Slider params values
        % slider.mu_sg   friction coefficient between slider and ground
        % slider.mu_sp   friction coefficient between slider and pusher
        % slider.ywidth  width of the slider along x-direction [m]
        % slider.xwidth  width of the slider along y-direction [m]
        % slider.area    slider area [m^2]
        % slider.m       slider mass [kg]

        % system dimensions
        nx = 5;     % number of state variables
        nu = 2;     % number of inputs

        % symbolic model
        sym_model = struct;
    end
    
    methods
        
        % Constructor
        function self = PusherSliderModel(name, slider_parameters)
            self@casadi.Callback();
            self.slider_params=slider_parameters;
            self.slider_params.f_max = self.slider_params.mu_sg*self.slider_params.m*helper.g;
            self.slider_params.tau_max = helper.tau_max_func(self.slider_params.mu_sg, self.slider_params.m, helper.g, self.slider_params.area, self.slider_params.xwidth, self.slider_params.ywidth);
            self.slider_params.c_ellipse = self.slider_params.tau_max/self.slider_params.f_max;
            construct(self, name);
        end

        % Contact surface
        function contact = contact_surface(self,S_p_x,S_p_y)
            %This function evaluate if the Pusher at (S_p_x,S_p_y) is in contact with the
            %Slider
            if abs(S_p_y) > self.slider_params.ywidth/2 || (S_p_x ~= -self.slider_params.xwidth/2)
               contact = false;
            else
               contact = true;
            end 
        end
        
        % Evaluate model (numerical)
        function x_dot = eval_model(self,x,u)
            % Pusher-Slider non linear model x_dot = f(x,u)
            % Input: x = [x_s, y_s, theta_s (rad), r_y],
            %        u = [u_n, u_t],
            % Output: x_dot

            % State variables
            theta_s = x(3);    % rotation angle of the slider frame w.r.t. the world frame
            S_p_x = x(4);      % x-coordinate of the pusher position w.r.t. the slider frame S
            S_p_y = x(5);      % y-coordinate of the pusher position w.r.t. the slider frame S
            u_n = u(1);        % normal pusher velocity w.r.t. the slider frame S
            u_t = u(2);        % tangential pusher velocity w.r.t. the slider frame S
            
            % 2D Rotation matrix about z-axis
            R_z = helper.my_rotz(theta_s); R_z = R_z(1:2,1:2);
            
            % Motion Cone
            gamma_l = (self.slider_params.mu_sp*self.slider_params.c_ellipse^2-S_p_x*S_p_y+self.slider_params.mu_sp*S_p_x^2)/(self.slider_params.c_ellipse^2+S_p_y^2-self.slider_params.mu_sp*S_p_x*S_p_y);
            gamma_r = (-self.slider_params.mu_sp*self.slider_params.c_ellipse^2-S_p_x*S_p_y-self.slider_params.mu_sp*S_p_x^2)/(self.slider_params.c_ellipse^2+S_p_y^2+self.slider_params.mu_sp*S_p_x*S_p_y);
            v_l = [1 gamma_l]'; % left edge of the motion cone
            v_r = [1 gamma_r]'; % right edge of the motion cone


            % Evaluating u_t/u_n
            u_fract = u_t/u_n;
            
            % Check boundary of the motion cone
            if not(contact_surface(self,S_p_x,S_p_y))
                mode="NC";
                disp("Not in contact!");
            else
                if (u_fract <= gamma_l) && (u_fract >= gamma_r)
                    mode = 'ST';
                    disp("Motion Cone CHECK: Sticking Mode")
                elseif u_fract > gamma_l
                    mode = 'SL';
                    disp("Motion Cone CHECK: Sliding Left Mode")
                else 
                    mode = 'SR';
                    disp("Motion Cone CHECK: Sliding Right Mode")
                end
            end

            % Model matrices
            factor_matrix = 1/(self.slider_params.c_ellipse^2+S_p_x^2+S_p_y^2);
            Q = [self.slider_params.c_ellipse^2+S_p_x^2 S_p_x*S_p_y; S_p_x*S_p_y self.slider_params.c_ellipse^2+S_p_y^2];
            
            switch mode
                case 'ST'
                    P = eye(2);
                    b = [-S_p_y S_p_x]';
                    %disp('sticking mode')
                case 'SL'                   % the pusher slides on the object surface
                    P = [v_l zeros(2,1)];
                    b = [-S_p_y+gamma_l*S_p_x 0]';
                    %disp('sliding left mode')
                case 'SR'                   % the pusher slides on the object surface
                    P = [v_r zeros(2,1)];
                    b = [-S_p_y+gamma_r*S_p_x 0]';
                    %disp('sliding right mode')
                otherwise                   % the pusher is not in contact with the slider
                    %disp('no contact')
                    x_dot = [0 0 0 u_n u_t]';
                    return;
            end
            c = eye(2)-factor_matrix*(Q*P+[-S_p_y; S_p_x]*b');
            F = [R_z*factor_matrix*Q*P; factor_matrix*b'; c];
            x_dot = F*[u_n u_t]';
        end
        
        % Symbolic model (symbolic)
        function model = symbolic_model(self)
            % This method returns the symbolic expression of the dynamical
            % model x_dot = f(x,u)
            import casadi.*
                      
            % named symbolic variables
            x = SX.sym('x');         % x-coordinate slider w.r.t. world frame [m]
            y = SX.sym('y');         % y-coordinate slider w.r.t. world frame [m]
            theta = SX.sym('theta'); % rotation angle of the slider frame w.r.t. the world frame [rad]
            S_p_x = SX.sym('S_p_x'); % x-coordinate of the pusher position w.r.t. the slider frame S [m]
            S_p_y = SX.sym('S_p_y'); % y-coordinate of the pusher position w.r.t. the slider frame S [m]
            
            u_n = SX.sym('u_n');         % normal pusher velocity w.r.t. slider frame S [m/s]
            u_t = SX.sym('u_t');         % tangential pusher velocity w.r.t. slider frame S [m/s]
            
            %% (unnamed) symbolic variables
            sym_x = vertcat(x,y,theta,S_p_x,S_p_y);
            sym_xdot = SX.sym('xdot', self.nx, 1);
            sym_u = vertcat(u_n,u_t);
            
            %% dynamics
            sin_theta = sin(theta);
            cos_theta = cos(theta);
            c_ellipse = self.slider_params.c_ellipse;
            mu_sp = self.slider_params.mu_sp;
            factor_matrix = 1/(c_ellipse^2+S_p_x^2+S_p_y^2);
            
            % motion cone
            gamma_l = (mu_sp*c_ellipse^2-S_p_x*S_p_y+mu_sp*S_p_x^2)/(c_ellipse^2+S_p_y^2-mu_sp*S_p_x*S_p_y);
            gamma_r = (-mu_sp*c_ellipse^2-S_p_x*S_p_y-mu_sp*S_p_x^2)/(c_ellipse^2+S_p_y^2+mu_sp*S_p_x*S_p_y);
            v_l = [1 gamma_l]';     % left edge of the motion cone
            v_r = [1 gamma_r]';     % right edge of the motion cone
            u_fract = u_t/u_n;
            
            % Model 
            
            R_z = [cos_theta -sin_theta;sin_theta cos_theta];
            Q = [c_ellipse^2+S_p_x^2 S_p_x*S_p_y; S_p_x*S_p_y c_ellipse^2+S_p_y^2];
            % Sticking
            P_st = eye(2);
            b_st = [-S_p_y S_p_x]';
            c_st = eye(2)-factor_matrix*(Q*P_st+[-S_p_y; S_p_x]*b_st');
            F_st = [R_z*factor_matrix*Q*P_st; factor_matrix*b_st'; c_st];
            x_dot_st = F_st*[u_n u_t]';
            % Sliding left
            P_sl = [v_l zeros(2,1)];
            b_sl = [-S_p_y+gamma_l*S_p_x 0]';
            c_sl = eye(2)-factor_matrix*(Q*P_sl+[-S_p_y; S_p_x]*b_sl');
            F_sl = [R_z*factor_matrix*Q*P_sl; factor_matrix*b_sl'; c_sl];
            x_dot_sl = F_sl*[u_n u_t]';
            % Sliding right
            P_sr = [v_r zeros(2,1)];
            b_sr = [-S_p_y+gamma_r*S_p_x 0]';
            c_sr = eye(2)-factor_matrix*(Q*P_sr+[-S_p_y; S_p_x]*b_sr');
            F_sr = [R_z*factor_matrix*Q*P_sr; factor_matrix*b_sr'; c_sr];
            x_dot_sr = F_sr*[u_n u_t]';
            
            
            % f = Function('f',{sym_x,sym_u},{(u_fract>=gamma_r)*x_dot_st*(u_fract<=gamma_l) ...
            %     + (u_fract>gamma_l)*x_dot_sl...
            %     + (u_fract<gamma_r)*x_dot_sr} ...
            %     );
            
            % explicit dynamic function
            expr_f_expl = vertcat((u_fract>=gamma_r)*x_dot_st*(u_fract<=gamma_l) ...
            + (u_fract>gamma_l)*x_dot_sl...
            + (u_fract<gamma_r)*x_dot_sr);
            
            % implicit dynamic function
            expr_f_impl = expr_f_expl - sym_xdot;
          
            %% populate structure
            model.nx = self.nx;
            model.nu = self.nu;
            model.sym_x = sym_x;
            model.sym_xdot = sym_xdot;
            model.sym_u = sym_u;
            model.expr_f_expl = expr_f_expl;
            model.expr_f_impl = expr_f_impl;  

            self.sym_model = model;
        end

    end
end