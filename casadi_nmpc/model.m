% Pusher-Slider non linear model x_dot = f(x,u)
% Input: x = [x_s, y_s, theta_s (rad), r_y], 
%        u = [u_n, u_t],
%        parameters = struct of the pusher_slider model parameters
% Output: x_dot

function x_dot = model(x, u, c)
x_dot = [1 1 1 1 1]';
disp(c);
%     u_n = u(1);
%     u_t = u(2);
%     % Model parameters
%     c_ellipse = slider.c_ellipse;
%     
%     % 2D Rotation matrix about z-axis
%     R_z = helper.my_rotz(theta_s); R_z = R_z(1:2,1:2);
%    
%     % Model matrices
%     factor_matrix = 1/(c_ellipse^2+S_p_x^2+S_p_y^2);
%     Q = [c_ellipse^2+S_p_x^2 S_p_x*S_p_y; S_p_x*S_p_y c_ellipse^2+S_p_y^2];
%     
%     % Motion Cone check
%     [mode, gamma_l, gamma_r] = motion_cone(u_n,u_t,x,slider);
%     
%     v_l = [1 gamma_l]'; % left edge of the motion cone
%     v_r = [1 gamma_r]'; % right edge of the motion cone
    
    
%     switch mode
%         case 'ST'                   
%             P = eye(2);
%             b = [-S_p_y S_p_x]';
% %             d = [0 0]';
% %             c = zeros(2,1);
%             disp('sticking mode')
%         case 'SL'                   % the pusher slides on the object surface
%             P = [v_l zeros(2,1)];
%             b = [-S_p_y+gamma_l*S_p_x 0]';
% %             d = zeros(2,1);
% %             c = [-gamma_l 0]';
%             disp('sliding left mode')
%         case 'SR'                   % the pusher slides on the object surface
%             P = [v_r zeros(2,1)];
%             b = [-S_p_y+gamma_r*S_p_x 0]';
% %             d = zeros(2,1);
% %             c = [-gamma_r 0]';
%             disp('sliding right mode')
%         otherwise                   % the pusher is not in contact with the slider
%             disp('no contact')
%             x_dot = [0 0 0 u_n u_t]';
%             return;
%     end
%     c = eye(2)-factor_matrix*(Q*P+[-S_p_y; S_p_x]*b');
%     F = [R_z*factor_matrix*Q*P; factor_matrix*b'; c];
%     x_dot = F*[u_n u_t]';

end
% import casadi.*
% classdef PusherSliderModel < casadi.Callback
%    
%   properties
%     slider_parameters
%     x_S
%     y_S
%   end
%   
%   methods
%     function self = PusherSliderModel(name, slider_parameters, x_S, y_S)
%       self@casadi.Callback();
%       self.slider_parameters = slider_parameters;
%       
%       % State variables
%       self.x_S = x_S;    % x-coordinate of the slider w.r.t. the world frame
%       self.y_S = y_S;    % y-coordinate of the slider w.r.t. the world frame
% 
%       construct(self, name);
%       x_S
%     end
% 
%     % Number of inputs and outputs
%     function v=get_n_in(self)
%       v=1;
%     end
%     function v=get_n_out(self)
%       v=1;
%     end
% 
%     % Initialize the object
%     function init(self)
%       disp('initializing object')
%     end
% 
%   end
% end