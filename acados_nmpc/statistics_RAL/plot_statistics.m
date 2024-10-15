load("NMPC/S-TRAJ/1_traj1_1_1_1.mat")
time = params.t;
traj = params.controller.y_ref;
end_array = length(traj);



%% Evaluate J
type = "NMPC";


W_x = 0.01*diag([100 100 0.1 0]);  % State matrix weight
W_x_e = 200*diag([1000 1000 0.1 0]); %diag([100 100 0 0 0]);
W_u = diag([1e-3 1e-3]);

t_dist = [10.0  550];

x0_x = [0.0260    0.0107    0.0155    0.0146   -0.0065];
x0_y = [0.0093   -0.0197    0.0124   -0.0281   -0.0134];
x0_theta = [-8.0492   -4.4300    0.9376    9.1501    9.2978]; %degree
x0_Spy = [-0.0382  -0.0248   -0.0010   -0.0070    0.0011];
y_dist = [0.0315    0.0406   -0.0273    0.0    0.0132];


numSim = length(y_dist);
filename_array = strings;
index_file_name = 1;
current_sim_num = 25;
% end_array = 200;
for index_t_dist = 2 : length(t_dist)
    for index_x0 = 1 : length(x0_x)
        for index_y_dist = 1:length(y_dist)
            current_sim_num = current_sim_num + 1;
            filename_array(index_file_name) = string(strcat(type,"/S-TRAJ/",num2str(current_sim_num),"_traj1_",num2str(index_t_dist),"_",num2str(index_x0),"_", num2str(index_y_dist)));
            index_file_name = index_file_name + 1;
        end
    end
end

% Evaluate index J for all data
num_sims = length(filename_array);
J = size(filename_array);
for i = 1:25%num_sims
    load(filename_array(i));
    if strcmp(type,"FOM")
        J(i) = evalJ([params.x_S(1:end_array)'; params.y_S(1:end_array)'; params.theta_S(1:end_array)';  [0;params.S_p_y(1:end_array-1)]']', traj(1:4,1:end_array)',W_x);
    else
        J(i) = evalJ([params.x_S(1:end_array)' params.y_S(1:end_array)' params.theta_S(1:end_array)' params.S_p_y(1:end_array)'], traj(1:4,:)',W_x);
    end
end

figure(1), hold on,
plot(J)
%% PLOT DEBUG
% type = "FOM";


load(strcat("FOM/STRAIGHT_LINE/1_traj1_1_1_1.mat"))
cad_model_path = "/home/workstation/pusher_slider_matlab/cad_models/cad_santal_centered_scaled_rotated_reduced.stl";
end_array = 199; %length(params.x_S)
helper.my_plot(time(1:end_array), [traj(1:3,1:end_array); traj(4,1:end_array)], params.x_S(1:end_array), params.y_S(1:end_array), params.theta_S(1:end_array), [params.S_p_y(1:end_array)], params.u_n(1:end_array), params.u_t(1:end_array));%, helper.convert_str2num(params.mode_vect));

load(strcat("NMPC/STRAIGHT_LINE/1_traj1_1_1_1.mat"))
helper.my_plot(time(1:end_array), [traj(1:3,1:end_array); traj(4,1:end_array)], params.x_S(1:end_array), params.y_S(1:end_array), params.theta_S(1:end_array), [params.S_p_y(1:end_array)], params.u_n(1:end_array), params.u_t(1:end_array));%, helper.convert_str2num(params.mode_vect));


% params.S_p_x = repmat(-0.034,1,length(params.x_S));
helper.my_animate(params.x_S,params.y_S,params.theta_S,params.S_p_x,params.S_p_y,time,0.1, [traj repmat(traj(:,end),1,(abs(length(params.x_S) - length(traj))))],cad_model_path);



function J = evalJ(x,xref,Wx)
    J = 0;
    for i = 1:length(x)
        J = J + (xref(i,:) - x(i,:))*Wx*(xref(i,:) - x(i,:))';
    end
end