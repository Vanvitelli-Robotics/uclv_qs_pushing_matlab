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


% %%%%%%%%%%%%%%%%%%%%%% SETUP REAL MODEL %%%%%%%%%%%%%%%%%%%%%%

% Pusher_Slider struct parameters
slider.mu_sg = 0.32;                                  % friction coefficient between slider and ground
slider.mu_sp = 0.19;                                 % friction coefficient between slider and pusher
slider.xwidth = 0.068;                               % width of the slider along x-direction [m]
slider.ywidth = 0.082;                                % width of the slider along y-direction [m]
slider.area = slider.xwidth * slider.ywidth;          % slider area [m^2]
slider.m = 0.2875;                                    % slider mass [kg]
plant_time_delay = 0;                               % delay of the plant [s]

% Create Pusher Slider object
cad_model_path = "../cad_models/cuboide_santal_resampled.stl";
order_spline = 3;
z_limit = 0.1;
p = PusherSliderModel('real_plant',slider, plant_time_delay,cad_model_path,order_spline,z_limit);
p.symbolic_model();


% %%%%%%%%%%%%%%%%%%%%%% SETUP NOMINAL PLANT %%%%%%%%%%%%%%%%%%%%%%

% Pusher_Slider struct parameters
slider_nominal.mu_sg = 0.32 - 0.1;                                  % friction coefficient between slider and ground
slider_nominal.mu_sp = 0.19 - 0.1;                                 % friction coefficient between slider and pusher
slider_nominal.xwidth = 0.068;                               % width of the slider along x-direction [m]
slider_nominal.ywidth = 0.082;                                % width of the slider along y-direction [m]
slider_nominal.area = slider.xwidth * slider.ywidth;          % slider area [m^2]
slider_nominal.m = 0.2875;                                    % slider mass [kg]
plant_time_delay_inc = 0.3;                               % delay of the plant [s]

% Create Pusher Slider object
% p_inc = PusherSliderModel('nominal_plant',slider_nominal, plant_time_delay_inc);
% p_inc.symbolic_model();

% %%%%%%%%%%%%%%%%%%%%%% SETUP CONTROLLER %%%%%%%%%%%%%%%%%%%%%%

% Controller parameters
sample_time = 0.05;
Hp = 25;

% Setup Controller and Optimization Object
controller = NMPC_controller('NMPC',p,sample_time,Hp);
controller.create_ocp_solver();


%% SETTING PARAMETERS FOR CONTROLLER AND PLANT

% Change delay of the plant and the delay to compensate with the controller
p.set_delay(0.35);
controller.set_delay_comp(0.35);

% Set initial condition
x0 = [0.0 0 deg2rad(0) slider.ywidth/2*0.3]';
controller.initial_condition_update(x0);

% Set matrix weights
% W_x = 10*diag([10 10 .001 0]);  % State matrix weight
% W_x_e = 100*diag([10 10 .001 0]); %diag([100 100 0 0 0]);
% W_u = diag([.1 1]);            % Control matrix weight
% 
W_x = 0.01*diag([100 100 .001 0]);  % State matrix weight
W_x_e = 200*diag([.05 .05 200 0]); %diag([100 100 0 0 0]);
W_u = diag([.1 0.8]);            % Control matrix weight


% controller.update_cost_function(W_x,W_u,W_x_e,Hp,Hp);
controller.update_cost_function(W_x,W_u,W_x_e,1,Hp-1);

% Set constraints
u_n_lb = 0.0; u_n_ub = 0.03;
u_t_lb = -0.03; u_t_ub = 0.03;
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
% x0_w = [0 0 0];
xf_w = [ 
      0.1 0.0 0;
        0.15 -0.1 0;
    0.25 -0.2 0;
    ];
traj_gen.waypoints_ = [x0_w; xf_w];
[time, traj] = traj_gen.waypoints_gen;
traj = [traj(1:3,:); traj(end,:)];
% traj = [traj repmat(traj(:,end),1,500)];
% time = [time(1:end-1) time(end):sample_time:(time(end)+sample_time*500)];
time_sim = time(end) + 2;

% Set control reference
u_n_ref = u_n_ub/2; u_t_ref = 0;
control_ref = repmat([u_n_ref; u_t_ref],1,length(time));

% Set overall reference
controller.set_reference_trajectory([traj; control_ref]);
% controller.set_reference_trajectory(traj);

% controller.create_ocp_solver();
% SIMULATION START

% If you want to simulate set simulation_ true and then set the
% type of simulation (simulink, matlab or real robot)
simulation_ = true;
sym_type = "matlab";
print_robot = false;


% Simulation
if simulation_ == true
    if(strcmp(sym_type,"simulink"))
        disp("SIMULINK SIMULATION")
        [x_s, y_s, theta_s, S_p_x, S_p_y, u_n, u_t, time_plot] = helper.closed_loop_simulink(time_sim);
        params = helper.save_parameters("exp1",[x_s; y_s; theta_s; S_p_x; S_p_y],[u_n; u_t],time_plot);

    elseif(strcmp(sym_type,"matlab"))
        disp("MATLAB SIMULATION")
        noise_ = true;
        debug_ = false;
        print_ = false;
        disturbance_ = true;
         
        [x_s, y_s, theta_s, S_p_x, S_p_y, u_n, u_t, time_plot,mode_vect, found_sol] = helper.closed_loop_matlab(p,controller,x0,time_sim,print_,noise_,debug_, disturbance_);
        params = helper.save_parameters("exp1_smooth_traj_",[x_s; y_s; theta_s; S_p_y],[u_n; u_t],time_plot, mode_vect);

    elseif(strcmp(sym_type,"robot"))
        disp("ROBOT EXPERIMENT")
        robot_main;
    else
        disp("Simulation type not valid!")
    end


end

%% PLOT
if sym_type == "robot"
    params.t = params.t(1:end-1);
    helper.my_plot_robot(params.t, [controller.y_ref(1:4,:) repmat(controller.y_ref(1:4,end),1,(abs(length(params.x_S) - length(controller.y_ref))))], params.x_S, params.y_S, params.theta_S, params.S_p_y, params.u_n, params.u_t);
else
    helper.my_plot(params.t, [controller.y_ref(1:3,:); controller.y_ref(4,:)], params.x_S, params.y_S, params.theta_S, params.S_p_y, params.u_n, params.u_t, controller.cost_function_vect, helper.convert_str2num(params.mode_vect));
end

%% ANIMATE
% helper.my_animate(params.x_S,params.y_S,params.theta_S,params.S_p_x,params.S_p_y, params.t,0.1, [controller.y_ref repmat(controller.y_ref(:,end),1,(abs(length(params.x_S) - length(controller.y_ref))))]);
params.S_p_x = repmat(-p.slider_params.xwidth/2,1,length(params.x_S));
helper.my_animate(params.x_S,params.y_S,params.theta_S,params.S_p_x,params.S_p_y, params.t,0.1, [controller.y_ref repmat(controller.y_ref(:,end),1,(abs(length(params.x_S) - length(controller.y_ref))))]);






