close all
clc

import casadi.*

% Let's define a function x_dot = f(x,u) with x = [x1 x2] and u = u1
x1 = MX.sym('x1');
x2 = MX.sym('x2');
x = [x1; x2];

u = MX.sym('u');

model = [(1-x2^2)*x1-x2+u; x1];
f = Function('f',{x,u},{model},{'x','u'},{'model'});

% Integration of function
intg_options = struct;
intg_option.tf = 1/2;
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

%MPC Parameters
Hp = 20; % prediction horizon
Hu = Hp; % control horizon
num_states = 2; % number of state variables
num_controls = 1; % number of control variables

% Optimal control problem
opti = casadi.Opti();
x = opti.variable(num_states, Hp+1); % state variables
u = opti.variable(num_controls, Hu);
p = opti.parameter(2,1); % initial condition

opti.minimize(sumsqr(x)+sumsqr(u));

for k=1:Hp
    opti.subject_to(x(:,k+1)==F(x(:,k),u(:,k)));
end

opti.subject_to(-1<=u<=1);
opti.subject_to(x(:,1)==p);


% Solver
% p_opts = struct('expand',true);
% s_opts = struct('max_iter',100);
% opti.solver('ipopt', p_opts, s_opts);
opti.solver('sqpmethod',struct('qpsol','qrqp'));
% opti.solver('ipopt');
% sol = opti.solve();

opti.set_value(p,[0;1]);
sol = opti.solve();

% PLOT
figure
hold on
tgrid = linspace(0,1/2,Hp+1);
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

x0 = [0;1];
x = x0;
for i=1:4*Hp
  u = full(M(x));

  U_log(:,i) = u;
  X_log(:,i) = x;

  % simulate system
  x = full(F(x,u)) + [0;rand*0.02];
end

