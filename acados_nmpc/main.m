clear all 
close all
clc

model_path = fullfile(pwd,'.');
addpath(model_path);
model_path = fullfile(pwd,'..');
addpath(model_path);
model_path = fullfile(pwd,'../cad_models');
addpath(model_path);

% Pusher_Slider struct parameters
slider.mu_sg = 0.32;                                  % friction coefficient between slider and ground
slider.mu_sp = 0.19;                                  % friction coefficient between slider and pusher
slider.ywidth = 0.0625;                               % width of the slider along x-direction [m]
slider.xwidth = 0.082;                                % width of the slider along y-direction [m]
slider.area = slider.xwidth * slider.ywidth;          % slider area [m^2]
slider.m = 0.2875;                                    % slider mass [kg]

%% plant
p = PusherSliderModel('Plant_real',slider);
p.symbolic_model();

% Initial condition
x0 = [0 0 deg2rad(0) -p.slider_params.xwidth/2 -p.slider_params.ywidth/2*0.5]';

%% controller
controller = NMPC_controller('NMPC',p,x0);
controller.create_ocp_solver();


%% use these commands to set new weigths and new initial conditions 
x0 = [0.15 0 deg2rad(0) -p.slider_params.xwidth/2 p.slider_params.ywidth/2*0]';
W_x = diag([1 1 1 1e-2 1e-2]);
W_u = diag([1 1]);

controller.initial_condition_update(x0);
controller.update_cost_function(W_x,W_u);


%% Simulation
time_sim = 10;
sim('simulation_model_closed_loop',time_sim);
x_s = out.signals.values(:,1);
y_s = out.signals.values(:,2);
theta_s = out.signals.values(:,3);
S_p_x = out.signals.values(:,4); 
S_p_y = out.signals.values(:,5);
u_n = u_control.signals.values(:,1);
u_t = u_control.signals.values(:,2);
my_animate(x_s,y_s,theta_s,S_p_x,S_p_y,controller.T/controller.N)






