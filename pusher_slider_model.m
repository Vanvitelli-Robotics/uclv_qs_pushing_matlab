% Pusher-Slider non linear model x_dot = f(x,u)
% Input: x = [x_s, y_s, theta_s (rad), r_y], 
%        u = [u_n, u_t],
%        parameters = struct of the pusher_slider model parameters
% Output: x_dot

function x_dot = pusher_slider_model(u_n, u_t, x, slider)
    
    % State variables
    theta_s = x(3);    % rotation angle of the slider frame w.r.t. the world frame
    S_p_x = x(4);      % x-coordinate of the pusher position w.r.t. the slider frame S
    S_p_y = x(5);      % y-coordinate of the pusher position w.r.t. the slider frame S
    
    % Model parameters
    c_ellipse = slider.c_ellipse;
    
    % 2D Rotation matrix about z-axis
    R_z = helper.my_rotz(theta_s); R_z = R_z(1:2,1:2);
    
    
    % Model matrices
    factor_matrix = 1/(c_ellipse^2+S_p_x^2+S_p_y^2);
    Q = factor_matrix*[c_ellipse^2+S_p_x^2 S_p_x*S_p_y; S_p_x*S_p_y c_ellipse^2+S_p_y^2];
    
    % Motion Cone check
    [mode, gamma_l, gamma_r] = motion_cone(u_n,u_t,x,slider);
    
    v_l = [1 gamma_l]'; % left edge of the motion cone
    v_r = [1 gamma_r]'; % right edge of the motion cone
    
    
    switch mode
        case 'ST'                   
            P = eye(2);
            b = [-S_p_y S_p_x]';
            d = [0 0]';
            c = zeros(2,1);
            disp('sticking mode')
        case 'SL'                   % the pusher slides on the object surface
            P = [v_l zeros(2,1)];
            b = [-S_p_y+gamma_l*S_p_x 0]';
            d = zeros(2,1);
            c = [-gamma_l 0]';
            disp('sliding left mode')
        case 'SR'                   % the pusher slides on the object surface
            P = [v_r zeros(2,1)];
            b = [-S_p_y+gamma_r*S_p_x 0]';
            d = zeros(2,1);
            c = [-gamma_r 0]';
            disp('sliding right mode')
        otherwise                   % the pusher is not in contact with the slider
            disp('no contact')      
    end
    
    F = [R_z*Q*P; factor_matrix*b'; d'; c'];
    x_dot = F*[u_n u_t]';

end