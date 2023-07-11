
classdef PusherSliderModel < casadi.Callback
    
    properties
        g = 9.81; % gravity acceleration [m/s^2]
        
        % Pusher_Slider parameters
        slider_params = struct;
        gamma_l;
        gamma_r;
    end
    
    methods
        
        % Constructor
        function self = PusherSliderModel(name, slider_parameters)
            self@casadi.Callback();
            self.slider_params=slider_parameters;
            self.slider_params.f_max = self.slider_params.mu_sg*self.slider_params.m*self.g;
            self.slider_params.tau_max = helper.tau_max_func(self.slider_params.mu_sg, self.slider_params.m, self.g, self.slider_params.area, self.slider_params.xwidth, self.slider_params.ywidth);
            self.slider_params.c_ellipse = self.slider_params.tau_max/self.slider_params.f_max;
            construct(self, name);
        end
        
        function x_dot = eval_model(self,mode,x,u)
            % Pusher-Slider non linear model x_dot = f(x,u)
            % Input: x = [x_s, y_s, theta_s (rad), r_y],
            %        u = [u_n, u_t],
            %        parameters = struct of the pusher_slider model parameters
            % Output: x_dot
            % State variables
            theta_s = x(3);    % rotation angle of the slider frame w.r.t. the world frame
            S_p_x = x(4);      % x-coordinate of the pusher position w.r.t. the slider frame S
            S_p_y = x(5);      % y-coordinate of the pusher position w.r.t. the slider frame S
            u_n = u(1);
            u_t = u(2);
            
            % 2D Rotation matrix about z-axis
            R_z = helper.my_rotz(theta_s); R_z = R_z(1:2,1:2);
            
            % Motion Cone
            self.gamma_l = (self.slider_params.mu_sp*self.slider_params.c_ellipse^2-S_p_x*S_p_y+self.slider_params.mu_sp*S_p_x^2)/(self.slider_params.c_ellipse^2+S_p_y^2-self.slider_params.mu_sp*S_p_x*S_p_y);
            self.gamma_r = (-self.slider_params.mu_sp*self.slider_params.c_ellipse^2-S_p_x*S_p_y-self.slider_params.mu_sp*S_p_x^2)/(self.slider_params.c_ellipse^2+S_p_y^2+self.slider_params.mu_sp*S_p_x*S_p_y);
            v_l = [1 self.gamma_l]'; % left edge of the motion cone
            v_r = [1 self.gamma_r]'; % right edge of the motion cone
            
            % Model matrices
            factor_matrix = 1/(self.slider_params.c_ellipse^2+S_p_x^2+S_p_y^2);
            Q = [self.slider_params.c_ellipse^2+S_p_x^2 S_p_x*S_p_y; S_p_x*S_p_y self.slider_params.c_ellipse^2+S_p_y^2];
            
            switch mode
                case 'ST'
                    P = eye(2);
                    b = [-S_p_y S_p_x]';
                    disp('sticking mode')
                case 'SL'                   % the pusher slides on the object surface
                    P = [v_l zeros(2,1)];
                    b = [-S_p_y+self.gamma_l*S_p_x 0]';
                    disp('sliding left mode')
                case 'SR'                   % the pusher slides on the object surface
                    P = [v_r zeros(2,1)];
                    b = [-S_p_y+self.gamma_r*S_p_x 0]';
                    disp('sliding right mode')
                otherwise                   % the pusher is not in contact with the slider
                    disp('no contact')
                    x_dot = [0 0 0 u_n u_t]';
                    return;
            end
            c = eye(2)-factor_matrix*(Q*P+[-S_p_y; S_p_x]*b');
            F = [R_z*factor_matrix*Q*P; factor_matrix*b'; c];
            x_dot = F*[u_n u_t]';
        end
    end
end