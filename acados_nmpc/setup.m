% This script realizes the setup for the nonlinear mpc controller with
% acados
% 

% Create Pusher Slider object
p = PusherSliderModel('pusher_slider_model',slider);
p.symbolic_model();

% Setup Controller and Optimization Object 
controller = NMPC_controller('NMPC',p,x0,linux_set);
controller.create_ocp_solver();
