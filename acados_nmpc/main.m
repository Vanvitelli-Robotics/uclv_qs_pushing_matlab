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

% Set initial condition
x0 = [0.15 0 deg2rad(0) -slider.xwidth/2 slider.ywidth/2*0]';

% Create optimization problem
setup

%% SIMULATION PART 
% If you want to simulate on simulink set simulation_ true
simulation_ = true;
time_sim = 10;

% Set initial condition
x0 = [0 0 deg2rad(0) -slider.xwidth/2 slider.ywidth/2*0.7]';

% Set matrix weights
W_x = diag([1 1 1 1e-2 1e-2]);  % State matrix weight
W_u = diag([1 1]);              % Control matrix weight

% Set initial condition and cost function 
controller.initial_condition_update(x0);
controller.update_cost_function(W_x,W_u);

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
    my_animate(x_s,y_s,theta_s,S_p_x,S_p_y,controller.T/controller.N)
end


    







