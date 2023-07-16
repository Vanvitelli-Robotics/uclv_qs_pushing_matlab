
time_exp = time_sim; % Time of the experiment [s]

rosinit('192.168.2.94',11310)

has_robot = false;
has_mpc_pose = false;

% tf object
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
x = x0;
u = [];
time_vec = [];
mpc_state_sub = rossubscriber("/mpc_pose","geometry_msgs/Pose2D", {@get_mpc_state,tftree,T_BS0, p, controller, x, u start_time, command_pub, command_msg, time_vec}, "BufferSize",1);


function get_mpc_state(mpc_state_sub, mpc_pose, tftree, T_BS0, p, controller, x, u, start_time, command_pub, command_msg, time_vec)
    %     tic
    t = rostime('now')-start_time;
    time_vec = [time_vec; t];

    if t.seconds > time_exp
        command_msg.Linear.X = 0;
        command_msg.Linear.Y = 0;
        send(command_pub,command_msg)
        pause(2);
        rosshutdown

        % Save experiment parameters
        params = helper.save_parameters("exp1",x,u,time_vec);
        return
    end

    % Get homogeneous transform from base_link to push_frame
    T_PB_ = getTransform(tftree, "base_link", "push_frame",rostime(0));
    T_PB = transformToMatrix(T_PB_);

    % Get homogeneous transform from push_frame to slider0
    T_PS0 = T_BS0 * T_PB;
    
    % Get the pusher position in the slider frame
    S0_p = [T_PS0(1,end);T_PS0(2,end)]; % position of the pusher w.r.t. the slider0 frame
    S0_s = [mpc_pose.X;mpc_pose.Y]; % position of the slider w.r.t. the slider0 frame
    S0_ps = S0_p-S0_s; % relative position between slider and pusher w.r.t. the slider0 frame
    S_p = helper.my_rotz_2d(mpc_pose.Theta)* S0_ps;

    % Get state [x y theta rx ry]
    x_ = [mpc_pose.X mpc_pose.Y mpc_pose.Theta S_p]';
    x = [x x_];

    % Delay Buffer Simulation
    xk = controller.delay_buffer_sim(p, x(:,end));
    u = [u controller.solve(xk)];

    controller.u_buff_contr = [u(:,end) controller.u_buff_contr(:,1:end-1)];

    % Send command to robot
    T_S0B = inv(T_BS0);
    u_robot = T_S0B(1:3,1:3)*[vipi;0];

    command_msg.Linear.X = u_robot(1);
    command_msg.Linear.Y = u_robot(2);
    send(command_pub,command_msg)
    %     toc
end

function T = transformToMatrix(transf)
    R = quat2rotm([transf.Transform.Rotation.W, transf.Transform.Rotation.X, transf.Transform.Rotation.Y, transf.Transform.Rotation.Z]);
    P = [transf.Transform.Translation.X, transf.Transform.Translation.Y, transf.Transform.Translation.Z]';

    T = [R P; 0 0 0 1];
end



