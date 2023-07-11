% Pusher-Slider non linear model x_dot = f(x,u)
% Input: x = [x_s, y_s, theta_s (rad), r_y], 
%        u = [u_n, u_t],
%        parameters = struct of the pusher_slider model parameters
% Output: x_dot

function x_dot = pusher_slider_model_sym(u_n, u_t, x, f)
   x_dot = full(f(x,[u_n u_t]))';
end