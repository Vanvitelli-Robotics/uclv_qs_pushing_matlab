clear all
close all
clc

% %%%%%%%%%%%%%%%%%%%%%% SETUP ACADOS %%%%%%%%%%%%%%%%%%%%%%
% Setting path variables
model_path = fullfile(pwd,'.');
addpath(model_path);
model_path = fullfile(pwd,'..');
addpath(model_path);
model_path = fullfile(pwd,'../cad_models');
addpath(model_path);

% Specify if linux or windows (true = linux, false = windows)
linux_set = false;

if linux_set == false
    env_vars_acados;
end


% %%%%%%%%%%%%%%%%%%%%%% SETUP MODEL %%%%%%%%%%%%%%%%%%%%%%

% Pusher_Slider struct parameters
slider.mu_sg = 0.32;                                  % friction coefficient between slider and ground
slider.mu_sp = 0.19;                                  % friction coefficient between slider and pusher
slider.ywidth = 0.0625;                               % width of the slider along x-direction [m]
slider.xwidth = 0.082;                                % width of the slider along y-direction [m]
slider.area = slider.xwidth * slider.ywidth;          % slider area [m^2]
slider.m = 0.2875;                                    % slider mass [kg]
plant_time_delay = 0.1;                               % delay of the plant [s]

% Create Pusher Slider object
p = PusherSliderModel('pusher_slider_model',slider, plant_time_delay);
p.symbolic_model();

% %%%%%%%%%%%%%%%%%%%%%% SETUP CONTROLLER %%%%%%%%%%%%%%%%%%%%%%

% Controller parameters
sample_time = 0.05;
Hp = 20;

% Setup Controller and Optimization Object
controller = NMPC_controller('NMPC',p,sample_time,Hp);
controller.create_ocp_solver();

%% SETTING PARAMETERS FOR CONTROLLER AND PLANT

% Change delay of the plant and the delay to compensate with the controller
p.set_delay(0.1);
controller.set_delay_comp(0.1);

% Set initial condition
x0 = [0 0 deg2rad(0) -slider.xwidth/2 slider.ywidth/2*0.7]';
controller.initial_condition_update(x0);

% Set matrix weights
W_x = diag([10 10 .1 0 0.0]);  % State matrix weight
W_x_e = 10*W_x;                %diag([100 20 .5 0 0]);
W_u = diag([1 10]);            % Control matrix weight
controller.update_cost_function(W_x,W_u,W_x_e,Hp,Hp);
controller.update_cost_function(W_x,W_u,W_x_e,1,Hp-1);

% Set constraints
u_n_lb = 0; u_n_ub = 0.05;
u_t_lb = -0.05; u_t_ub = 0.05;
controller.update_constraints(u_n_ub, u_t_ub, u_n_lb, u_t_lb);

% Create desired trajectory
xf = [0.3 0.03 0 x0(4) 0]';
xf(3) = acos((xf(1)-x0(1))/(norm(xf(1:2)-x0(1:2))));

traj_gen = TrajectoryGenerator(sample_time,u_n_ub/2);
traj_gen.set_plot = false;

time_sim = 20;
t0 = 0; tf = time_sim*1;
traj_gen.set_target(x0,xf,t0,tf);
% [time, traj] = traj_gen.straight_line(true);

x0_w = [x0(1:2)' 0];
xf_w = [0.3 0.01 0;
    0.3 0.02 0
    ];
traj_gen.waypoints_ = [x0_w; xf_w];
[time, traj] = traj_gen.waypoints_gen;

% Set control reference
u_n_ref = 0; u_t_ref = 0;
control_ref = [u_n_ref; u_t_ref].*zeros(controller.sym_model.nu,length(time));

% Set overall reference
controller.set_reference_trajectory([traj; control_ref]);


%% SIMULATION START

% If you want to simulate on simulink set simulation_ true and then set the
% type of simulation (simulink, matlab or real robot)
simulation_ = true;
sym_type = "matlab";


% Simulation
if simulation_ == true
    if(strcmp(sym_type,"simulink"))
        [x_s, y_s, theta_s, S_p_x, S_p_y, u_n, u_t, time_plot] = helper.closed_loop_simulink(time_sim);
        params = helper.save_parameters("exp1",[x_s; y_s; theta_s; S_p_x; S_p_y],[u_n; u_t],time_plot);

    elseif(strcmp(sym_type,"matlab"))
        [x_s, y_s, theta_s, S_p_x, S_p_y, u_n, u_t, time_plot] = helper.closed_loop_matlab(p,controller,x0,time_sim,true);
        params = helper.save_parameters("exp1",[x_s; y_s; theta_s; S_p_x; S_p_y],[u_n; u_t],time_plot);

    elseif(strcmp(sym_type,"robot"))
        robot_main;
    else
        disp("Simulation type not valid!")
    end
    helper.my_plot(params.t, traj, params.x_S, params.y_S, params.theta_S, params.S_p_x, params.S_p_y, params.u_n, params.u_t)
    helper.my_animate(params.x_S,params.y_S,params.theta_S,params.S_p_x,params.S_p_y,controller.sample_time, traj)

    %     helper.my_plot(time_plot, traj, x_s, y_s, theta_s, S_p_x, S_p_y, u_n, u_t)
    %     helper.my_animate(x_s,y_s,theta_s,S_p_x,S_p_y,controller.sample_time, traj)

end









