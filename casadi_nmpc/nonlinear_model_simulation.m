clear variables
close all
clc

% Helpfull variables
g = 9.81; % gravity acceleration [m/s^2]

% Pusher_Slider struct parameters
slider.mu_sg = 0.32;                                  % friction coefficient between slider and ground
slider.mu_sp = 0.19;                                  % friction coefficient between slider and pusher
slider.ywidth = 0.0625;                               % width of the slider along x-direction [m]
slider.xwidth = 0.082;                                % width of the slider along y-direction [m]
slider.area = slider.xwidth * slider.ywidth;          % slider area [m^2]
slider.m = 0.2875;                                    % slider mass [kg]
% slider.f_max = slider.mu_sg*slider.m*g;               % maximum force with ellipsoidal approximation of LS [N]
% slider.tau_max = helper.tau_max_func(slider.mu_sg, slider.m, g, slider.area, slider.xwidth, slider.ywidth); % maximum torque with ellipsoidal approximation of LS [Nm]
% slider.c_ellipse = slider.tau_max/slider.f_max;

% Initial condition
x0 = [0 0 deg2rad(0) -slider.xwidth/2 -0.02]';

%% Simulation
time_sim = 100;
sim('simulation_model',time_sim);
x_s = out.signals.values(:,1);
y_s = out.signals.values(:,2);
theta_s = out.signals.values(:,3);
S_p_x = out.signals.values(:,4); 
S_p_y = out.signals.values(:,5);
my_animate(x_s,y_s,theta_s,S_p_x,S_p_y)

