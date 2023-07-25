function x_dot = eval_model(x,u,params)
    % This method can be used to evaluate the pusher-Slider non linear model x_dot = f(x,u)
    % Input: x = [x_s, y_s, theta_s (rad), r_y],
    %        u = [u_n, u_t],
    % Output: x_dot

    % State variables
    theta_s = x(3);    % rotation angle of the slider frame w.r.t. the world frame
    S_p_x = -0.034;% x(4);      % x-coordinate of the pusher position w.r.t. the slider frame S
    S_p_y = x(4);      % y-coordinate of the pusher position w.r.t. the slider frame S
    u_n = u(1);        % normal pusher velocity w.r.t. the slider frame S
    u_t = u(2);        % tangential pusher velocity w.r.t. the slider frame S

    % 2D Rotation matrix about z-axis
    R_z = helper.my_rotz(theta_s); R_z = R_z(1:2,1:2);

    % Motion Cone edges
    gamma_l = (params.mu_sp*params.c_ellipse^2-S_p_x*S_p_y+params.mu_sp*S_p_x^2)/(params.c_ellipse^2+S_p_y^2-params.mu_sp*S_p_x*S_p_y);
    gamma_r = (-params.mu_sp*params.c_ellipse^2-S_p_x*S_p_y-params.mu_sp*S_p_x^2)/(params.c_ellipse^2+S_p_y^2+params.mu_sp*S_p_x*S_p_y);
    v_l = [1 gamma_l]'; % left edge of the motion cone
    v_r = [1 gamma_r]'; % right edge of the motion cone


    % Evaluating u_t/u_n
    u_fract = u_t/u_n;

    % Check boundary of the motion cone
%     if not(contact_surface(self,S_p_x,S_p_y))
%         mode="NC";
%         disp("Not in contact!");
%     else
        if (u_fract <= gamma_l) && (u_fract >= gamma_r)
            mode = "ST";
            %                     disp("Motion Cone CHECK: Sticking Mode")
        elseif u_fract > gamma_l
            mode = "SL";
            %                     disp("Motion Cone CHECK: Sliding Left Mode")
        else
            mode = "SR";
            %                     disp("Motion Cone CHECK: Sliding Right Mode")
        end
%     end


    % Model matrices
    factor_matrix = 1/(params.c_ellipse^2+S_p_x^2+S_p_y^2);
    Q = [params.c_ellipse^2+S_p_x^2 S_p_x*S_p_y; S_p_x*S_p_y params.c_ellipse^2+S_p_y^2];
    %             mode = 'ST';
    switch mode
        case 'ST'
            P = eye(2);
            b = [-S_p_y S_p_x]';
            %                     disp('sticking mode')
        case 'SL'                   % the pusher slides on the object surface
            P = [v_l zeros(2,1)];
            b = [-S_p_y+gamma_l*S_p_x 0]';
            %                     disp('sliding left mode')
        case 'SR'                   % the pusher slides on the object surface
            P = [v_r zeros(2,1)];
            b = [-S_p_y+gamma_r*S_p_x 0]';
            %                     disp('sliding right mode')
        otherwise                   % the pusher is not in contact with the slider
            %disp('no contact')
            x_dot = [0 0  u_n u_t]';
            return;
    end
    c = eye(2)-factor_matrix*(Q*P+[-S_p_y; S_p_x]*b');
    F = [R_z*factor_matrix*Q*P; factor_matrix*b'; c(end,:)];
    x_dot = F*[u_n u_t]';
end