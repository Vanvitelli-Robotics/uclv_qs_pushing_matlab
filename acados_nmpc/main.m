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
x0 = [0 0 deg2rad(0) -p.slider_params.xwidth/2 p.slider_params.ywidth/2*0.9]';

% u = [0.0 -0.0]';
% 
% p.eval_model(x0,u);

%% controller
y_ref = [0.05 0 deg2rad(0) -p.slider_params.xwidth/2 0]';

controller = NMPC_controller('NMPC',p);
controller.create_ocp_solver();



%controller.solve(x0,y_ref);


%% Simulation
time_sim = 30;
sim('simulation_model_closed_loop',time_sim);
x_s = out.signals.values(:,1);
y_s = out.signals.values(:,2);
theta_s = out.signals.values(:,3);
S_p_x = out.signals.values(:,4); 
S_p_y = out.signals.values(:,5);

my_animate(x_s,y_s,theta_s,S_p_x,S_p_y)





