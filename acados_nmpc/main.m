clear all
close all
clc

% Setting path variables
model_path = fullfile(pwd,'.');
addpath(model_path);
model_path = fullfile(pwd,'..');
addpath(model_path);
model_path = fullfile(pwd,'../cad_models');
addpath(model_path);

% Specify if linux or windows (true = linux, false = windows)
linux_set = true;

% Pusher_Slider struct parameters
slider.mu_sg = 0.32;                                  % friction coefficient between slider and ground
slider.mu_sp = 0.19;                                  % friction coefficient between slider and pusher
slider.ywidth = 0.0625;                               % width of the slider along x-direction [m]
slider.xwidth = 0.082;                                % width of the slider along y-direction [m]
slider.area = slider.xwidth * slider.ywidth;          % slider area [m^2]
slider.m = 0.2875;                                    % slider mass [kg]

% Create Pusher Slider object
p = PusherSliderModel('pusher_slider_model',slider);
p.symbolic_model();

% Controller parameters
sample_time = 0.05;
Hp = 20;
u_n_lb = 0; u_n_ub = 0.5;
u_t_lb = -0.5; u_t_ub = 0.5;
y_lb = -p.slider_params.ywidth/2;
y_ub = p.slider_params.ywidth/2;

% Setup Controller and Optimization Object 
controller = NMPC_controller('NMPC',p,linux_set,sample_time,Hp,u_n_ub, u_t_ub, u_n_lb, u_t_lb);
controller.create_ocp_solver();

%% SIMULATION PART 
close all

% If you want to simulate on simulink set simulation_ true
simulation_ = true;
time_sim = 10;

% Set initial condition
x0 = [0 0 deg2rad(0) -slider.xwidth/2 slider.ywidth/2*0.7]';
controller.initial_condition_update(x0);

% Set matrix weights
W_x = diag([1 1 .01 0 0]);  % State matrix weight
W_x_e = 10*W_x;%diag([10 20 .5 0 0]);
W_u = diag([1 1]);              % Control matrix weight
controller.update_cost_function(W_x,W_u,W_x_e);

% Set constraints
u_n_lb = 0; u_n_ub = 0.05;
u_t_lb = -0.05; u_t_ub = 0.05;
controller.update_constraints(u_n_ub, u_t_ub, u_n_lb, u_t_lb);

% Create desired trajectory
xf = [0.3 0.06 0 x0(4) 0]'; 
xf(3) = acos((xf(1)-x0(1))/(norm(xf(1:2)-x0(1:2))));
t0 = 0; tf = time_sim;
traj_gen = TrajectoryGenerator(sample_time);
traj_gen.set_plot = false;
traj_gen.set_target(x0,xf,t0,tf);
[time, traj] = traj_gen.straight_line;

% Set control reference
u_n_ref = 0; u_t_ref = 0;
control_ref = [u_n_ref; u_t_ref].*zeros(controller.sym_model.nu,length(time));

% Set overall reference
controller.set_reference_trajectory([traj; control_ref]);

% Simulation 
if simulation_ == true
    sim('simulation_model_closed_loop',time_sim);
    x_s = out.signals.values(:,1);
    y_s = out.signals.values(:,2);
    theta_s = out.signals.values(:,3);
    S_p_x = out.signals.values(:,4);
    S_p_y = out.signals.values(:,5);
    u_n = u_control.signals.values(:,1);
    u_t = u_control.signals.values(:,2);
    helper.my_animate(x_s,y_s,theta_s,S_p_x,S_p_y,controller.sample_time, traj)
    helper.my_plot(out.time, traj, x_s, y_s, theta_s, S_p_x, S_p_y, u_n, u_t)
end


    







