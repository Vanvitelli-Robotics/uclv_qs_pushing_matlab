% Motion Cone constraints check
% Input: x = [x_s, y_s, theta_s (rad), r_y], 
%        u = [u_n, u_t],
%        parameters = struct of the pusher_slider model parameters
% Output: mode, gamma_l, gamma_r

function [mode, gamma_l, gamma_r] = motion_cone(u_n,u_t,x,slider)
    S_p_x = x(4);      % x-coordinate of the pusher position w.r.t. the slider frame S
    S_p_y = x(5);      % y-coordinate of the pusher position w.r.t. the slider frame S
    mu_sp = slider.mu_sp;
    c_ellipse = slider.c_ellipse;
    
    % Calculating edges of motion cone
%     gamma_l = (mu_sp*c_ellipse^2+S_p_x*S_p_y-mu_sp*S_p_x^2)/(c_ellipse^2-S_p_y^2+mu_sp*S_p_x*S_p_y);
%     gamma_r = (-mu_sp*c_ellipse^2+S_p_x*S_p_y+mu_sp*S_p_x^2)/(c_ellipse^2-S_p_y^2-mu_sp*S_p_x*S_p_y);

    gamma_l = (mu_sp*c_ellipse^2-S_p_x*S_p_y+mu_sp*S_p_x^2)/(c_ellipse^2+S_p_y^2-mu_sp*S_p_x*S_p_y);
    gamma_r = (-mu_sp*c_ellipse^2-S_p_x*S_p_y-mu_sp*S_p_x^2)/(c_ellipse^2+S_p_y^2+mu_sp*S_p_x*S_p_y);
    
    
    % Evaluating u_t/u_n
    u_fract = u_t/u_n;
    
    % Check boundary of the motion cone
    if abs(S_p_y) > slider.ywidth/2 || (S_p_x ~= -slider.xwidth/2)
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
end
