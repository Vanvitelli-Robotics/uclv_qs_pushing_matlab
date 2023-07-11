% This script solve the non linear problem 
clear all
close all
clc

% acados initialitation
env_vars_acados
check_acados_requirements()

model_path = fullfile(pwd,'..');
addpath(model_path);
%% model dynamics
model = pusher_slider_model_MPC;
nx = model.nx;
nu = model.nu;
ny = size(model.cost_expr_y, 1);      % used in simulink example
ny_e = size(model.cost_expr_y_e, 1);

%% discretization
N = 20;
sample_time = 0.050; %controller sample time [s]
T = N*sample_time;   % time horizon length [s]

% initial condition
x0 = [0; 0; deg2rad(0); -model.slider.xwidth/2; 0.02];

nlp_solver = 'sqp';                      % sqp, sqp_rti
qp_solver = 'partial_condensing_hpipm';  % full_condensing_hpipm, partial_condensing_hpipm, full_condensing_qpoases, full_condensing_daqp
qp_solver_cond_N = 5;                    % for partial condensing
% integrator type
sim_method = 'erk';                      % erk, irk, irk_gnsf

%% model to create the solver
ocp_model = acados_ocp_model();
model_name = 'pusher_slider';

%% acados ocp model
ocp_model.set('name', model_name);
ocp_model.set('T', T);

% symbolics
ocp_model.set('sym_x', model.sym_x);
ocp_model.set('sym_u', model.sym_u);
ocp_model.set('sym_xdot', model.sym_xdot);

% cost
ocp_model.set('cost_expr_ext_cost', model.expr_ext_cost);
ocp_model.set('cost_expr_ext_cost_e', model.expr_ext_cost_e);

% dynamics
if (strcmp(sim_method, 'erk'))
    ocp_model.set('dyn_type', 'explicit');
    ocp_model.set('dyn_expr_f', model.expr_f_expl);
else % irk irk_gnsf
    ocp_model.set('dyn_type', 'implicit');
    ocp_model.set('dyn_expr_f', model.expr_f_impl);
end

% constraints
ocp_model.set('constr_type', 'auto');
ocp_model.set('constr_expr_h', model.expr_h);

% on the variables defined to be constrained
h_max = [model.slider.ywidth/2 0.05 0.05];
h_min = [-model.slider.ywidth/2 0 -0.05];
ocp_model.set('constr_lh', U_min); % lower bound on h
ocp_model.set('constr_uh', U_max);  % upper bound on h

% % on the state
% ocp_model.set('constr_Jbx',eye(nx));
% ocp_model.set('constr_lbx',[-100 -100 -100 -100 -model.slider.ywidth/2]);
% ocp_model.set('constr_ubx',[100 100 100 100 model.slider.ywidth/2]);
% 
% % inputs 
% ocp_model.set('constr_Jbu',eye(nu));
% ocp_model.set('constr_lbu',U_min);
% ocp_model.set('constr_ubu',U_max);

ocp_model.set('constr_x0', x0);
% ... see ocp_model.model_struct to see what other fields can be set

%% acados ocp set opts
ocp_opts = acados_ocp_opts();
ocp_opts.set('param_scheme_N', N);
ocp_opts.set('nlp_solver', nlp_solver);
ocp_opts.set('sim_method', sim_method);
ocp_opts.set('qp_solver', qp_solver);
ocp_opts.set('qp_solver_cond_N', qp_solver_cond_N);
ocp_opts.set('ext_fun_compile_flags', '-O2'); % '-O2'
ocp_opts.set('compile_model','false');
%ocp_opts.set('codgen_model','false');
ocp_opts.set('compile_interface','false');
ocp_opts.set('output_dir',fullfile(pwd,'build'));
% ... see ocp_opts.opts_struct to see what other fields can be set

%% create ocp solver
ocp = acados_ocp(ocp_model, ocp_opts);

x_traj_init = zeros(nx, N+1);
u_traj_init = zeros(nu, N);

%% call ocp solver
% update initial state
ocp.set('constr_x0', x0);

% reference
ocp.set('cost_y_ref_e', [0.05 0.05 rad2deg(0) -model.slider.xwidth/2 0]);


% set trajectory initialization
ocp.set('init_x', x_traj_init);
ocp.set('init_u', u_traj_init);
ocp.set('init_pi', ones(nx, N))

% change values for specific shooting node using:
%   ocp.set('field', value, optional: stage_index)
ocp.set('constr_lbx', x0, 0)

% solve
tic
ocp.solve();
toc
% get solution
utraj = ocp.get('u');
xtraj = ocp.get('x');

status = ocp.get('status'); % 0 - success
ocp.print('stat')

%% Plots
ts = linspace(0, T, N+1);
figure; hold on;
States = {'x', 'y', 'theta', 'S_p_x', 'S_p_y'};
for i=1:length(States)
    subplot(length(States), 1, i);
    plot(ts, xtraj(i,:)); grid on;
    ylabel(States{i});
    xlabel('t [s]')
end

figure
Inputs = {'u_n', 'u_t'};
for i=1:length(Inputs)
    subplot(length(States), 1, i);
    plot(ts, [utraj(i,:),utraj(i,end)]); grid on;
    ylabel(Inputs{i});
    xlabel('t [s]')
end
grid on

%% go embedded
% to generate templated C code
% ocp.generate_c_code;
