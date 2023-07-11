% This function returns the model of the pusher-slider system used by the
% NMPC. 
% The returned struct "model" contains:
%       - system dimensions
%       - system parameters
%       - symbolic variables
%       - dynamics model expression (symbolic)
%       - constraints of the plant
%       - cost matrices

function model = pusher_slider_model_MPC()

import casadi.*

%% system dimensions
nx = 5;
nu = 2;

%% system parameters
g = 9.81;                                             % gravity acceleration [m/s^2]
slider.mu_sg = 0.32;                                  % friction coefficient between slider and ground
slider.mu_sp = 0.19;                                  % friction coefficient between slider and pusher
slider.ywidth = 0.0625;                               % width of the slider along x-direction [m]
slider.xwidth = 0.082;                                % width of the slider along y-direction [m]
slider.area = slider.xwidth * slider.ywidth;          % slider area [m^2]
slider.m = 0.2875;                                    % slider mass [kg]
slider.f_max = slider.mu_sg*slider.m*g;               % maximum force with ellipsoidal approximation of LS [N]
slider.tau_max = helper.tau_max_func(slider.mu_sg, slider.m, g, slider.area, slider.xwidth, slider.ywidth); % maximum torque with ellipsoidal approximation of LS [Nm]
slider.c_ellipse = slider.tau_max/slider.f_max;

%% named symbolic variables
x = SX.sym('x');         % x-coordinate slider w.r.t. world frame [m]
y = SX.sym('y');         % y-coordinate slider w.r.t. world frame [m]
theta = SX.sym('theta'); % rotation angle of the slider frame w.r.t. the world frame [rad]
S_p_x = SX.sym('S_p_x'); % x-coordinate of the pusher position w.r.t. the slider frame S [m]
S_p_y = SX.sym('S_p_y'); % y-coordinate of the pusher position w.r.t. the slider frame S [m]

u_n = SX.sym('u_n');         % normal pusher velocity w.r.t. slider frame S [m/s]
u_t = SX.sym('u_t');         % tangential pusher velocity w.r.t. slider frame S [m/s]

%% (unnamed) symbolic variables
sym_x = vertcat(x,y,theta,S_p_x,S_p_y);
sym_xdot = SX.sym('xdot', nx, 1);
sym_u = vertcat(u_n,u_t);

%% dynamics
sin_theta = sin(theta);
cos_theta = cos(theta);
c_ellipse = slider.c_ellipse;
mu_sp = slider.mu_sp;
factor_matrix = 1/(c_ellipse^2+S_p_x^2+S_p_y^2);

% motion cone
gamma_l = (mu_sp*c_ellipse^2-S_p_x*S_p_y+mu_sp*S_p_x^2)/(c_ellipse^2+S_p_y^2-mu_sp*S_p_x*S_p_y);
gamma_r = (-mu_sp*c_ellipse^2-S_p_x*S_p_y-mu_sp*S_p_x^2)/(c_ellipse^2+S_p_y^2+mu_sp*S_p_x*S_p_y);
v_l = [1 gamma_l]';     % left edge of the motion cone
v_r = [1 gamma_r]';     % right edge of the motion cone
u_fract = u_t/u_n;

% model 

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


%% cost
W_x = diag([1, 1, 1, 1, 1]);
W_u = diag([1 1]);
expr_ext_cost_e = sym_x'* W_x * sym_x;
expr_ext_cost = expr_ext_cost_e + sym_u' * W_u * sym_u;
% nonlinear least sqares
cost_expr_y = vertcat(sym_x, sym_u);
W = blkdiag(W_x, W_u);
model.cost_expr_y_e = sym_x;
model.W_e = W_x;

%% constraints
expr_h = vertcat(S_p_y,sym_u);
%% populate structure
model.nx = nx;
model.nu = nu;
model.sym_x = sym_x;
model.sym_xdot = sym_xdot;
model.sym_u = sym_u;
model.expr_f_expl = expr_f_expl;
model.expr_f_impl = expr_f_impl;
model.expr_h = expr_h;
model.expr_ext_cost = expr_ext_cost;
model.expr_ext_cost_e = expr_ext_cost_e;

model.cost_expr_y = cost_expr_y;
model.W = W;
model.slider = slider;

end
