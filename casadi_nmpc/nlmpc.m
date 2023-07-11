close all
clc
% addpath("/home/workstation/casadi-3.6.3-linux64-matlab2018b");
addpath("C:\Users\saraf\casadi-3.6.3-windows64-matlab2018b");

%%
import casadi.*

% State variables
x_S = MX.sym('x_S');    % x-coordinate of the slider w.r.t. the world frame
y_S = MX.sym('y_S');    % y-coordinate of the slider w.r.t. the world frame
theta_S = MX.sym('theta_S');    % rotation angle of the slider frame w.r.t. the world frame
S_p_x = MX.sym('S_p_x');    % x-coordinate of the pusher position w.r.t. the slider frame S
S_p_y = MX.sym('S_p_y');    % y-coordinate of the pusher position w.r.t. the slider frame S
x = [x_S y_S theta_S S_p_x S_p_y]';

% Control variables
u_n = MX.sym('u_n');
u_t = MX.sym('u_t');
u = [u_n u_t]';
u_fract = u_t/u_n;

obj = PusherSliderModel('obj', slider);

x_dot_ST = obj.eval_model('ST',x,u);
x_dot_SR = obj.eval_model('SR',x,u);
x_dot_SL = obj.eval_model('SL',x,u);

f = Function('f',{x,u},{(u_fract<=obj.gamma_l)*x_dot_ST*(u_fract>=obj.gamma_r)...
    + (u_fract>obj.gamma_l)*x_dot_SL...
    + (u_fract<obj.gamma_r)*x_dot_SR},{'x','u'},{'x_dot'});


%%
%MPC Parameters
Hp = 20; % prediction horizon
Hu = Hp; % control horizon
num_states = 5; % number of state variables
num_controls = 2; % number of control variables
T = 0.05; %sample time

% Integration of function
intg_options = struct;
intg_option.tf = Hp*T;
intg_option.simplify = true;
intg_options.number_of_finite_elements = 4;

dae = struct;
dae.x = x; % state
dae.p = u;
dae.ode = f(x,u);

intg = integrator('intg','rk',dae,intg_options);
res = intg('x0',x,'p',u);
x_new = res.xf;

F = Function('F',{x,u},{x_new},{'x','u'},{'x_next'});
%%
% Optimal control problem
opti = casadi.Opti();
x = opti.variable(num_states, Hp+1); % state variables
u = opti.variable(num_controls, Hu);
p = opti.parameter(5,1); % initial condition

opti.minimize(sumsqr(x)+sumsqr(u));

for k=1:Hp
    opti.subject_to(x(:,k+1)==F(x(:,k),u(:,k)));
end




% Solver
% p_opts = struct('expand',true);
% s_opts = struct('max_iter',100);
% opti.solver('ipopt', p_opts, s_opts);
opti.solver('sqpmethod',struct('qpsol','qrqp'));
% opti.solver('ipopt');
% sol = opti.solve();

opti.set_value(p,x0);
sol = opti.solve();
%%
% PLOT
figure
hold on
tgrid = linspace(0,0.05,Hp+1);
plot(tgrid,sol.value(x));
stairs(tgrid, [sol.value(u) nan], '-.');
xlabel('t [s]');
legend('x1','x2','u');

%% MPC Loop


% Make the solver silent
opts = struct;
opts.qpsol = 'qrqp';
opts.print_header = false;
opts.print_iteration = false;
opts.print_time = false;
opts.qpsol_options.print_iter = false;
opts.qpsol_options.print_header = false;
opts.qpsol_options.print_info = false;
opti.solver('sqpmethod',opts);



M = opti.to_function('M',{p},{u(:,1)},{'p'},{'u_opt'});

% MPC loop
X_log = [];
U_log = [];

x = x0;
for i=1:4*Hp
  u = full(M(x));

  U_log(:,i) = u;
  X_log(:,i) = x;

  % simulate system
  x = full(F(x,u)) + [0;rand*0.02;rand*0.02;rand*0.02;rand*0.02];
end

