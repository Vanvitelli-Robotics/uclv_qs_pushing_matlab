
time_exp = time_sim; % Time of the experiment [s]


rosshutdown
rosinit('192.168.2.94',11310)
% ros.Node("/robot_test",'192.168.2.94',11310);

has_robot = false;
has_mpc_pose = false;

% tf object
clear tftree
tftree = rostf;
tftree.BufferTime = 0.5;
pause(2);

% Get homogeneous transform from base_link to slider0
T_BS0_ = getTransform(tftree, "slider0", "base_link",rostime(0));
T_BS0 = transformToMatrix(T_BS0_);

command_pub = rospublisher("/command_vel_des","geometry_msgs/Twist");
command_msg = rosmessage(command_pub);
start_time = rostime('now');

% State and control variables
global x x_sim u time_vec print_ found_sol f s0_spline options

x = [];
x_sim = [];
u = [];
time_vec = [];
found_sol = [];
print_ = print_robot;
params = struct;

f = @(x,c) norm([c.S_p_x c.S_p_y] - c.plant.SP.evalSpline(c.plant.SP.FC,x))^2;
options = optimoptions('fminunc','Display','off');
s0_spline = 0;


mpc_state_sub = rossubscriber("/mpc_state","std_msgs/Float64MultiArray", {@get_mpc_state,T_BS0, tftree, p, controller, start_time, command_pub, command_msg, time_exp, params}, "BufferSize",1);


function get_mpc_state(mpc_state_sub, mpc_state, T_BS0, tftree, p, controller, start_time, command_pub, command_msg, time_exp, params)
    global x x_sim u time_vec print_ found_sol f s0_spline options
    tic
    mpc_pose.X = mpc_state.Data(1);
    mpc_pose.Y = mpc_state.Data(2);
    mpc_pose.Theta = mpc_state.Data(3);
    t = rostime('now')-start_time;
    time_vec = [time_vec; t.Sec + 10^-9*t.Nsec];
    disp(t.Sec + 10^-9*t.Nsec)

    if t.seconds > time_exp

        command_msg.Linear.X = 0;
        command_msg.Linear.Y = 0;
        send(command_pub,command_msg)
        delete(mpc_state_sub);
        delete(command_pub);
        clear tftree
%         rosshutdown
        disp("Saving parameters")
        helper.save_parameters("exp_robot_traj",x,x_sim,u,time_vec-time_vec(1), params);
        return
    end


    % Get the pusher position in the slider frame
    S0_s = [mpc_pose.X;mpc_pose.Y]; % position of the pusher w.r.t. the slider0 frame
    S0_p = [mpc_state.Data(4);mpc_state.Data(5)]; % position of the slider w.r.t. the slider0 frame
    S0_ps = S0_p-S0_s; % relative position between slider and pusher w.r.t. the slider0 frame
    S_p = helper.my_rotz_2d(-mpc_pose.Theta)* S0_ps;

    % Get the curvilinear abscissa from S_p by solving an optimization problem
    c.S_p_x = S_p(1); c.S_p_y = S_p(2); c.plant = p;
    smin_spline = fminunc(@(x) f(x,c),s0_spline,options);
    s0_spline = smin_spline;

   
    % Get state [x y theta rx ry]
%     x_ = [mpc_pose.X mpc_pose.Y mpc_pose.Theta S_p(2)]';
%     x = [x x_];

    % Get state [x y theta s]
    x_ = [mpc_pose.X mpc_pose.Y mpc_pose.Theta smin_spline]';
    x = [x x_];

    % Delay Buffer Simulation
    xk = controller.delay_buffer_sim(p, x(:,end));
    x_sim = [x_sim xk];
    u = [u controller.solve(xk,round(time_vec(end)/controller.sample_time)+controller.delay_buff_comp)];
%     if t.seconds > 1
%         u_ = [0.0; 0.005];
%     else
%         u_ = [0 0]';
%     end
%     u = [u u_];

    controller.u_buff_contr = [u(:,end) controller.u_buff_contr(:,1:end-1)];
    status = controller.ocp_solver.get('status');
    
    if print_ == true
        sqp_iter = controller.ocp_solver.get('sqp_iter');
        time_tot = controller.ocp_solver.get('time_tot');
        time_lin = controller.ocp_solver.get('time_lin');
        time_qp_sol = controller.ocp_solver.get('time_qp_sol');

        fprintf('\nstatus = %d, sqp_iter = %d, time_int = %f [ms] (time_lin = %f [ms], time_qp_sol = %f [ms])\n',...
            status, sqp_iter, time_tot*1e3, time_lin*1e3, time_qp_sol*1e3);
    end
    if status~=0
        disp('acados ocp solver failed');
        %             keyboard
        found_sol = [found_sol false];
    else
        found_sol = [found_sol true];
    end


    % Send command to robot
    T_S0B = inv(T_BS0);
    R_NT_S_ = p.SP.R_NT_fun(1);
    R_NT_S = full(R_NT_S_(mod(smin_spline,p.SP.b)));
    R_NT_S_3d = [R_NT_S zeros(2,1); zeros(1,2), 1];
    R_NT_W = T_S0B(1:3,1:3)*helper.my_rotz(xk(3))*R_NT_S_3d;
    u_robot = R_NT_W*[u(:,end);0];
    
    % Calculate omega-y to rotate the end-effector
%     R_P_W = getTransform(tftree, "base_link", "push_frame",rostime(0));


    command_msg.Linear.X = u_robot(1);
    command_msg.Linear.Y = u_robot(2);
    % command_msg.Angular.Y = 
    send(command_pub,command_msg)

    %     helper.save_parameters("exp_robot",x,u,time_vec, params);
    toc
end

function T = transformToMatrix(transf)
    R = quat2rotm([transf.Transform.Rotation.W, transf.Transform.Rotation.X, transf.Transform.Rotation.Y, transf.Transform.Rotation.Z]);
    P = [transf.Transform.Translation.X, transf.Transform.Translation.Y, transf.Transform.Translation.Z]';

    T = [R P; 0 0 0 1];
end




