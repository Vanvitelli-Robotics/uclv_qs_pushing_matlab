clear all
% close all
clc

% %%%%%%%%%%%%%%%%%%%%%% SETUP ACADOS %%%%%%%%%%%%%%%%%%%%%%
% Setting path variables
model_path = fullfile(pwd,'.');
addpath(model_path);
model_path = fullfile(pwd,'..');
addpath(model_path);
model_path = fullfile(pwd,'../cad_models');
addpath(model_path);
model_path = fullfile(pwd,'./objects_database');
addpath(model_path);

% Specify if linux or windows (true = linux, false = windows)
linux_set = true;

if linux_set == false
    env_vars_acados;
end


% %%%%%%%%%%%%%%%%%%%%%% SETUP REAL MODEL %%%%%%%%%%%%%%%%%%%%%%

slider = object_selection('santal');
plant_time_delay = 0;                               % delay of the plant [s]

% Create Pusher Slider object
cad_model_path = slider.cad_model_path; %"/home/workstation/pusher_slider_matlab/cad_models/cad_santal_centered_scaled_rotated_reduced.stl"; %"../cad_models/cuboide_santal_rotated.stl";
pcl_path = slider.pcl_path; %'/home/workstation/pusher_slider_matlab/cad_models/santal_planar_surface_simplified.ply';
order_spline = 3;
p = PusherSliderModel('real_plant',slider, plant_time_delay,cad_model_path,order_spline,pcl_path);
% p.symbolic_model();
p.symbolic_model_variable_shape();



% %%%%%%%%%%%%%%%%%%%%%% SETUP CONTROLLER %%%%%%%%%%%%%%%%%%%%%%

% Controller parameters
sample_time = 0.07;
Hp = 10;

% Setup Controller and Optimization Object
controller = NMPC_controller('NMPC',p,sample_time,Hp);
controller.create_ocp_solver();


%% SETTING PARAMETERS FOR CONTROLLER AND PLANT

% Change delay of the plant and the delay to compensate with the controller
p.set_delay(0.33); %0.35
controller.set_delay_comp(0.33);

% Set initial condition
% x0 = [0.0 0 deg2rad(0) slider.ywidth/2*0.3]';
x0 = [0 0 deg2rad(0) -0.05]';
controller.initial_condition_update(x0);

% Set matrix weights
% Working matrix 
% % W_x = 0.01*diag([100 100 .001 0]);  % State matrix weight
% % W_x_e = 200*diag([.05 .05 200 0]); %diag([100 100 0 0 0]);
% % W_u = diag([.1 0.8]);            % Control matrix weight

% Matrix for variable shape 
% % W_x = 0.01*diag([100 100 .001 0]);  % State matrix weight
% % W_x_e = 200*diag([1000 1000 100 0]); %diag([100 100 0 0 0]);
% % W_u = diag([0 0]);  

W_x = 0.01*diag([100 100 .001 0]);  % State matrix weight
W_x_e = 200*diag([1000 1000 10 0]); %diag([100 100 0 0 0]);
W_u = diag([0 0]);  

% controller.update_cost_function(W_x,W_u,W_x_e,Hp,Hp);
controller.update_cost_function(W_x,W_u,W_x_e,1,Hp-1);

% Set constraints
u_n_lb = 0.0; u_n_ub = 0.02;
u_t_lb = -0.03; u_t_ub = 0.03;
controller.update_constraints(u_n_ub, u_t_ub, u_n_lb, u_t_lb);

% set constraints tangential velocity 
controller.set_v_alpha(0.008*200);

% Create desired trajectory
xf = [0.3 0.03 0 x0(4) 0.07]';
xf(3) = acos((xf(1)-x0(1))/(norm(xf(1:2)-x0(1:2))));

traj_gen = TrajectoryGenerator(sample_time,u_n_ub/2);
traj_gen.set_plot = false;

time_sim = 10;
t0 = 0; tf = time_sim*1;
traj_gen.set_target(x0,xf,t0,tf);
% [time, traj] = traj_gen.straight_line(true);

x0_w = [x0(1:2)' 0];
% x0_w = [0 0 0];
xf_w = [ 
      0.2 0.0 0;
%         0.15 -0.1 0;
%     0.25 -0.2 0;
    ];
traj_gen.waypoints_ = [x0_w; xf_w];
[time, traj] = traj_gen.waypoints_gen;
traj = [traj(1:3,:); traj(end,:)];

time_sim = time(end) + 15;

% Set control reference
u_n_ref = u_n_ub/2; u_t_ref = 0;
control_ref = repmat([u_n_ref; u_t_ref],1,length(time));

% Set overall reference
controller.set_reference_trajectory([traj; control_ref]);
% controller.set_reference_trajectory(traj);


%%%%%%%%%%%%% SIMULATION START %%%%%%%%%%%%%%%%%%%%%
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
        disturbance_ = false;
         
        [x_s, y_s, theta_s, S_p_x, S_p_y, u_n, u_t, time_plot,mode_vect, found_sol] = helper.closed_loop_matlab(p,controller,x0,time_sim,print_,noise_,debug_, disturbance_);
%         [x_s, y_s, theta_s, S_p_x, S_p_y, u_n, u_t, time_plot] = helper.open_loop_matlab(p,x0,0.0,0.05, time_sim,sample_time, noise_);
        
        params = helper.save_parameters("exp1_open_loop",[x_s; y_s; theta_s; S_p_y],[u_n; u_t],time_plot);
        params.S_p_x = S_p_x;
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
    helper.my_plot(params.t, [controller.y_ref(1:3,:); controller.y_ref(4,:)], params.x_S, params.y_S, params.theta_S, params.S_p_y, params.u_n, params.u_t, controller.cost_function_vect);%, helper.convert_str2num(params.mode_vect));
end

%% ANIMATE
% helper.my_animate(params.x_S,params.y_S,params.theta_S,params.S_p_x,params.S_p_y, params.t,0.1, [controller.y_ref repmat(controller.y_ref(:,end),1,(abs(length(params.x_S) - length(controller.y_ref))))]);
% params.S_p_x = repmat(-p.slider_params.xwidth/2,1,length(params.x_S));
if sym_type == "robot"
    S_p = p.SP.evalSpline(p.SP.FC,params.S_p_y);
    params.S_p_x = S_p(:,1);
    params.S_p_y = S_p(:,2);
end
helper.my_animate(params.x_S,params.y_S,params.theta_S,params.S_p_x,params.S_p_y, params.t,0.1, [controller.y_ref repmat(controller.y_ref(:,end),1,(abs(length(params.x_S) - length(controller.y_ref))))],cad_model_path);






