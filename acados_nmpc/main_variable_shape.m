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
slider.ywidth = 0.068;                               % width of the slider along x-direction [m]
slider.xwidth = 0.085;                                % width of the slider along y-direction [m]
slider.area = slider.xwidth * slider.ywidth;          % slider area [m^2]
slider.m = 0.2875;                                    % slider mass [kg]
plant_time_delay = 0;                               % delay of the plant [s]

% Create Pusher Slider object
%cad_santal_centered1
%cuboide_santal
cad_model_path = "../cad_models/cuboide_santal_resampled.stl";
order_spline = 3;
z_limit = 0.1;
p = PusherSliderModel('pusher_slider_model',slider, plant_time_delay,cad_model_path,order_spline,z_limit);
p.symbolic_model_variable_shape();

%%
time_sim = 20;
sample_time = 0.05;
% x0 = [0 0 deg2rad(0) -slider.ywidth/2*0.8]';
x0 = [0 0 deg2rad(0) p.SP.b+0.01]';
% velocità pusher in terna normale/tangenziale 
u_n_ = 0.001;
u_t_ = 0.05;

[x_s, y_s, theta_s, S_p_x, S_p_y, u_n, u_t, time_plot] = helper.open_loop_matlab(p,x0,u_n_,u_t_,time_sim,sample_time,false);
% params = helper.save_parameters("exp1_no_noise",[x_s; y_s; theta_s; S_p_x; S_p_y],[u_n; u_t],time_plot);

traj = zeros(p.nx+p.nu,length(time_plot));
helper.my_animate(x_s, y_s, theta_s, S_p_x, S_p_y,time_plot,sample_time,traj)
