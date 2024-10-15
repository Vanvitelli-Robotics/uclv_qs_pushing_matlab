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
model_path = fullfile(pwd,'./objects_database');
addpath(model_path);

% Specify if linux or windows (true = linux, false = windows)
linux_set = true;

if linux_set == false
    env_vars_acados;
end

% %%%%%%%%%%%%%%%%%%%%%% SETUP REAL MODEL %%%%%%%%%%%%%%%%%%%%%%
% Select the object from object_database file
object_name = "santal";
slider = object_selection(object_name); % load the physical parameters of the selected object
plant_time_delay = 0;                   % delay of the plant [s]

% Create Pusher Slider object
cad_model_path = slider.cad_model_path;
pcl_path = slider.pcl_path; 
order_spline = 3;
p = PusherSliderModel('real_plant',slider, plant_time_delay,cad_model_path,order_spline,pcl_path,object_name);
p.symbolic_model_variable_shape();


% %%%%%%%%%%%%%%%%%%%%%% SETUP CONTROLLER %%%%%%%%%%%%%%%%%%%%%%
% Controller parameters
sample_time = 0.05;
Hp = 10;

% Setup Controller and Optimization Object
controller = NMPC_controller('NMPC',p,sample_time,Hp);
controller.create_ocp_solver();


%% SETTING PARAMETERS FOR CONTROLLER AND PLANT

% t_dist = [4/0.05 7/0.05];
t_dist = [15/0.05 27/0.05];

x0_x = [0.0260   0.0107    0.0155    0.0146   -0.0065];
x0_y = [0.0093   -0.0197    0.0124   -0.0081   -0.0134];
x0_theta = [-8.0492   -4.4300    0.9376    9.1501    9.2978]; %degree
x0_Spy = [-0.0382  -0.0248   -0.0010   -0.0070    0.0011];
y_dist = [0.0315    0.0406   -0.015   0.0001    0.0132];



numSim = length(y_dist);
current_sim_num = 0;

for index_t_dist = 1 : length(t_dist)
    if index_t_dist == 1
        y_dist = [-0.035   -0.0161   -0.0127   -0.0021   -0.0018];
    else
        y_dist = [0.0215    -0.0306   -0.015   0.017    0.0132];
    end

    for index_x0 = 1 : length(x0_x)
        for index_y_dist = 1: length(y_dist)
            % Change delay of the plant and the delay to compensate with the controller
            p.set_delay(0.35*0);
            controller.set_delay_comp(0.35*0);

            % Set initial condition
            x0 = [x0_x(index_x0) x0_y(index_x0) deg2rad(x0_theta(index_x0)) x0_Spy(index_x0)]';
            controller.initial_condition_update(x0);

            % Set matrix weights
            W_x = 0.01*diag([100 100 0.1 0]);  
            W_x_e = 200*diag([1000 1000 0.1 0]); 
            W_u = diag([1e-3 1e-3]);

            controller.update_cost_function(W_x,W_u,W_x_e,0,Hp-1);

            % Set constraints
            u_n_lb = 0.0; u_n_ub = 0.03;
            u_t_lb = -0.05; u_t_ub = 0.05;
            % controller.update_constraints(u_n_ub, u_t_ub, u_n_lb, u_t_lb);

            % set constraints tangential velocity
            % alpha santal = 0.005*200;
            % alpha montana = 0.0027*200;
            % controller.set_v_alpha(0.005*200);%0.1650);%0.035

            % Create desired trajectory
            xf = [0.3 0.03 0 x0(4) 0.07]';
            xf(3) = acos((xf(1)-x0(1))/(norm(xf(1:2)-x0(1:2))));

            traj_gen = TrajectoryGenerator(sample_time,0.01);
            traj_gen.set_plot = true;

            time_sim = 10;
            t0 = 0; tf = time_sim*1;
            traj_gen.set_target(x0,xf,t0,tf);
            % [time, traj] = traj_gen.straight_line(true);

            x0_w = [0 0 0];
            % x0_w = [0 0 0];
            % xf_w = [
            % %       0.05 0 0; %trajectory exp
            %         0.1 0.0 0;
            %         0.15 0.1 0;
            %         0.25 0.2 0;
            %     ];

            % xf_w = [
            % %       0.05 0 0;   %exp cluttered scene
            %         0.0 0.11 0;
            %         -0.05 0.113 0;
            %     ];
            %
            % xf_w = [
            % %       0.05 0 0;   %exp cluttered scene
            %         0.05 0.0 0;
            %         0.13 0 0;
            %     ];

            % disturbance trajectory
            %
            % xf_w = [
            % %       0.05 0 0; %trajectory exp
            %         0.20 0.0 0;
            %         0.20 0.05 0;
            %         0.20 0.15 0;
            %         0.20 0.20 0;
            %         0 0.20 0;
            %     ];

            % xf_w = [
            % %       0.05 0 0; %trajectory exp
            %         0.12 0.0 0;
            %         0.12 0.06 0;
            %         0.12 0.12 0;
            % %         0.20 0.20 0;
            % %         0 0.20 0;
            %     ];

% % % %             xf_w = [
% % % %                 %       0.05 0 0; %trajectory exp
% % % % %                 0.0 0.1 0;
% % % % %                 0.0 0.30 0;
% % % %                 %         0.20 0.20 0;
% % % %                 %         0 0.20 0;
% % % %                 0.10 0 0;
% % % %                 ];



% % % %             traj_gen.waypoints_ = [x0_w; xf_w];
            % traj_gen.waypoints_velocities = [0.005 0.01 0.01];
            % traj_gen.waypoints_velocities = [0.005 0.003];
% % % %             traj_gen.waypoints_velocities = [0.010];% 0.01];

            % disturbance parameters
            % traj_gen.waypoints_velocities = [0.005 0.003 0.005];% 0.003 0.005];

% % % % % %             [time, traj] = traj_gen.waypoints_gen;
%             [time, traj] = traj_gen.waypoint_gen_fixed_angle;
% % % % % %             traj = [traj(1:3,:); traj(end,:)];
            load('x_finals.mat');
            traj = [x_finals_struct.x; x_finals_struct.y; x_finals_struct.theta; zeros(1,length(x_finals_struct.x))];
            time_sim = length(x_finals_struct.x)*sample_time;
            time = 0:sample_time:time_sim;
            time = time(1:end-1);

            % time_sim = x_finals_struct.t(end);%time(end)+15;
            % time = x_finals_struct.t(1:length(x_finals_struct.x));

            % % load('time_demo_pushing_UTRAJ.mat','time')
            % % load('traj_demo_pushing_UTRAJ.mat','traj')
            % % traj(2,:) = -traj(2,:);

            %
            % load('time_demo_pushing.mat','time')
            % load('traj_demo_pushing.mat','traj')
            % time_sim = time(end)+10;


            % Set control reference
            u_n_ref = 0.0; u_t_ref = 0;
            control_ref = repmat([u_n_ref; u_t_ref],1,length(time));

            % Set overall reference
            controller.set_reference_trajectory([traj; control_ref]);


            %%%%%%%%%%%%% SIMULATION START %%%%%%%%%%%%%%%%%%%%%
            % If you want to simulate set simulation_ true and then set the
            % type of simulation (simulink, matlab or robot)
            simulation_ = true;
            sym_type = "matlab";
            print_robot = false;

            current_sim_num = current_sim_num +1;
            total_sim_num = length(t_dist)*length(y_dist)*length(x0_x);
            disp(strcat("START SIMULATION ",num2str(current_sim_num),"/",num2str(total_sim_num)))
            
            % Simulation
            if simulation_ == true
                if(strcmp(sym_type,"simulink"))
                    disp("SIMULINK SIMULATION")
                    [x_s, y_s, theta_s, S_p_x, S_p_y, u_n, u_t, time_plot] = helper.closed_loop_simulink(time_sim);
                    params = helper.save_parameters("exp1",[x_s; y_s; theta_s; S_p_x; S_p_y],[u_n; u_t],time_plot);

                elseif(strcmp(sym_type,"matlab"))
                    disp("MATLAB SIMULATION")
                    noise_ = false;
                    debug_ = false;
                    print_ = false;
                    disturbance_ = true;

                    [x_s, x_sim, y_s, theta_s, S_p_x, S_p_y, u_n, u_t, time_plot,mode_vect, found_sol] = helper.closed_loop_matlab(p,controller,x0,time_sim,print_,noise_,debug_, disturbance_,y_dist(index_y_dist),t_dist(index_t_dist));
                    %         [x_s, y_s, theta_s, S_p_x, S_p_y, u_n, u_t, time_plot] = helper.open_loop_matlab(p,x0,0.0,-0.05, time_sim,sample_time, noise_);
                    filename = strcat("/home/workstation/pusher_slider_matlab/acados_nmpc/statistics_RAL/NMPC/",num2str(current_sim_num),"_traj1_",num2str(index_t_dist),"_",num2str(index_x0),"_",num2str(index_y_dist));
                    params = helper.save_parameters(filename,[x_s; y_s; theta_s; S_p_y],x_sim, [u_n; u_t],time_plot);
                    params.S_p_x = S_p_x;
                elseif(strcmp(sym_type,"robot"))
                    disp("ROBOT EXPERIMENT")
                    robot_main;
                else
                    disp("Simulation type not valid!")
                end


            end




            % PLOT
            params.controller.y_ref = controller.y_ref(:,1:min(length(params.x_S),length(controller.y_ref(1,:))));
            params.controller.cost_function_vect = controller.cost_function_vect;
            params.S_p_x = repmat(-p.slider_params.xwidth/2,1,length(params.x_S));
            save(filename,"params")

            if sym_type == "robot"
                params.t = params.t(1:end-1);
                time_plot(1:length(controller.y_ref)),
                helper.my_plot_robot(params.t(1:end-1), time, controller.y_ref(1:4,1:length(time)), params.x_S, params.y_S, params.theta_S, params.S_p_y, params.u_n, params.u_t);
            else
                helper.my_plot(params.t(1:end), [params.controller.y_ref(1:3,:); params.controller.y_ref(4,:)], params.x_S, params.y_S, params.theta_S, params.S_p_y, params.u_n, params.u_t, params.controller.cost_function_vect);%, helper.convert_str2num(params.mode_vect));
            end
%                pause;
        end
    end
end


return

%% ANIMATE

if sym_type == "robot"
    S_p = p.SP.evalSpline(p.SP.FC,params.S_p_y);
    params.S_p_x = S_p(:,1);
    params.S_p_y = S_p(:,2);
end

helper.my_animate(params.x_S,params.y_S,params.theta_S,params.S_p_x,params.S_p_y, params.t,0.5, [controller.y_ref repmat(controller.y_ref(:,end),1,(abs(length(params.x_S) - length(controller.y_ref))))],cad_model_path);






