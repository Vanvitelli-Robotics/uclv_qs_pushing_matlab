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

% %%%%%%%%%%%%%%%%%%%%%% SETUP MODEL %%%%%%%%%%%%%%%%%%%%%%
g = 9.81;
% Pusher_Slider struct parameters
slider.mu_sg = 0.32;                                  % friction coefficient between slider and ground
slider.mu_sp = 0.19;                                 % friction coefficient between slider and pusher
slider.xwidth = 0.068;                               % width of the slider along x-direction [m]
slider.ywidth = 0.082;                                % width of the slider along y-direction [m]
slider.area = slider.xwidth * slider.ywidth;          % slider area [m^2]
slider.m = 0.2875;                                    % slider mass [kg]
plant_time_delay = 0;                               % delay of the plant [s]

% Create Pusher Slider object
p = PusherSliderModel('pusher_slider_model',slider, plant_time_delay);
slider.f_max = slider.mu_sg*slider.m*g;               % maximum force with ellipsoidal approximation of LS [N]
slider.tau_max = p.tau_max_func(slider.mu_sg, slider.m, g, slider.area, slider.xwidth, slider.ywidth); % maximum torque with ellipsoidal approximation of LS [Nm]
slider.c_ellipse = slider.tau_max/slider.f_max;
p.symbolic_model();

% %%%%%%%%%%%%%%%%%%%%%%%%% NONLINEAR MPC %%%%%%%%%%%%%%%%%%%%
nx = 4;
nu = 2;
ny = 4;

nlobj = nlmpc(nx,ny,nu);

Ts = 0.03;
Hp = 20;
Hu = 10;
nlobj.Ts = Ts;
nlobj.PredictionHorizon = Hp;
nlobj.ControlHorizon = Hu;

nlobj.Model.StateFcn = @eval_model_continuous; %"eval_model_continuous";
nlobj.Model.IsContinuousTime = true;
nlobj.Model.NumberOfParameters = 1;

% State setting
nlobj.States(1).Name = "x_s";
nlobj.States(1).Units = "m";

nlobj.States(2).Name = "y_s";
nlobj.States(2).Units = "m";

nlobj.States(3).Name = "theta_s";
nlobj.States(3).Units = "rad";

nlobj.States(4).Name = "S_p_y";
nlobj.States(4).Units = "m";
nlobj.States(4).Min = -0.9*p.slider_params.ywidth/2;
nlobj.States(4).Max = 0.9*p.slider_params.ywidth/2;


% Control setting
u_n_lb = 0.0; u_n_ub = 0.03;
u_t_lb = -0.05; u_t_ub = 0.05;

nlobj.ManipulatedVariables(1).Name = "u_n";
nlobj.ManipulatedVariables(1).Units = "m/s";
nlobj.ManipulatedVariables(1).Max = u_n_ub;
nlobj.ManipulatedVariables(1).Min = u_n_lb;
% nlobj.ManipulatedVariables(1).RateMax = 0.01;
% nlobj.ManipulatedVariables(1).RateMin = -0.01;

nlobj.ManipulatedVariables(2).Name = "u_t";
nlobj.ManipulatedVariables(2).Units = "m/s";
nlobj.ManipulatedVariables(2).Max = u_t_ub;
nlobj.ManipulatedVariables(2).Min = u_t_lb;
% nlobj.ManipulatedVariables(2).RateMax = 0.02;
% nlobj.ManipulatedVariables(2).RateMin = -0.02;


% Cost function
W_x = diag([1000 1000 .001 0]);
W_u = diag([0 0]);
W_deltau = [0.0001 0.0001];
nlobj.Weights.OutputVariables = W_x;
nlobj.Weights.ManipulatedVariables = W_u;
nlobj.Weights.ManipulatedVariablesRate = W_deltau;

% Solver Options
nlobj.Optimization.SolverOptions.ConstraintTolerance = 1e-10;
nlobj.Optimization.SolverOptions.StepTolerance = 1e-10;
% nlobj.Optimization.SolverOptions.Algorithm = "interior-point";
% nlobj.Optimization.SolverOptions.HonorBounds = true;
% nlobj.Optimization.SolverOptions.Display = "iter";
nlobj.Optimization.SolverOptions.EnableFeasibilityMode = true;
% nlobj.Optimization.SolverOptions.PlotFcn = {@optimplotx,@optimplotconstrviolation,@optimplotfval,@optimplotfirstorderopt};
nlobj.Optimization.SolverOptions.UseParallel = true;

%% TRAJECTORY
% Create desired trajectory
x0 = [0.0 0 deg2rad(0) slider.ywidth/2*0]';
xf = [0.3 0.03 0 x0(4) 0]';
xf(3) = acos((xf(1)-x0(1))/(norm(xf(1:2)-x0(1:2))));

traj_gen = TrajectoryGenerator(Ts,u_n_ub/2);
traj_gen.set_plot = false;
traj_gen.x0 = x0;


x0_w = [x0(1:2)' 0];
% x0_w = [0 0 0];
xf_w = [
    0.1 0.0 0;
%         0.15 0.1 0;
%         0.25 0.2 0;
    ];
traj_gen.waypoints_ = [x0_w; xf_w];
[time, traj] = traj_gen.waypoints_gen;
traj = [traj(1:3,:); traj(end,:)];

time_sim = time(end);


%% SIMULATION
use_mex_ = false;

% Time of overall simulation
time_sim_vec = 0:Ts:time_sim;
time_sim_ = length(time_sim_vec);

% State and control variables
x = zeros(nx, time_sim_+1);
x(:,1) = x0;
u = zeros(nu, time_sim_);
mode_vect = string(zeros(time_sim_,1));
found_sol = true(length(time_sim_vec),1);
cost = zeros(1, time_sim_);
cost_dbg = cost;

nloptions = nlmpcmoveopt;
nloptions.Parameters = {slider};

noise_ = false;
debug_ = false;
print_ = true;
disturbance_ = true;
u0 = [0 0]';
u_last = u0;
MVopt = zeros(nlobj.PredictionHorizon,nu);


if use_mex_ == true
    [coreData,onlineData] = getCodeGenerationData(nlobj,x0,u_last,nloptions.Parameters);
    mexFcn = buildMEX(nlobj,"myController",coreData,onlineData);
end

%             tic;
for i = 1:time_sim_
    disp("Time: " + i*Ts);
    if disturbance_ == true
        if i == 33
            disp("Disturbance")
            x(2,i) = x(2,i) + 0.01;
            x(4,i) = x(4,i) - 0.01;
        end
    end


    % noise simulation
    if(noise_ == true)
        x(:,i) = x(:,i) + [1e-5*randn(1,2) randn()*1e-3 0 randn()*1e-4]';
    end

    % solve NMPC
    nloptions = nlmpc_updateX0(nlobj,nloptions,x(:,i), MVopt);
    if i+nlobj.PredictionHorizon > size(traj,2)
        traj_ = [traj(:,i:end) repmat(traj(:,end),1,i+nlobj.PredictionHorizon-1 - size(traj,2))]';
    else
        traj_ = traj(:,i:i+nlobj.PredictionHorizon-1)';
    end
    if use_mex_ == true
        onlineData.X0 = nloptions.X0;
        onlineData.ref = traj_;
        tic;
        [u(:,i),onlineData, info] = myController(x(:,i),u_last,onlineData);
        toc;
    else
        tic;
        [u(:,i),nloptions,info] = nlmpcmove(nlobj,x(:,i),u_last,traj_,[],nloptions);
        toc;
    end


    cost_dbg(i) = evalFcst_nlmpc(nlobj,x(:,i), info.MVopt, u_last, traj(:,i)',[],nloptions);

    if false && (  i <= 3 || (i>=30 && i<=inf*35) )
        disp(['contour ' num2str(i)])
        un = u_n_lb:0.0005:u_n_ub;
        ut = u_t_lb:0.0005:u_t_ub;

        %         un = [u_n_lb:0.001:0.02 0.02:0.0005:u_n_ub];
        %         ut = [u_t_lb:0.0005:-0.01 -0.01:0.001:u_t_ub];

        [UN,UT] = meshgrid(un,ut);
        FF = zeros(size(UN));
        MV__ = info.MVopt;

        TRAJ__ = traj(:,i)';
        X__ = x(:,i);
        MV__2end = MV__(2:end,:);
        parfor jj=1:numel(UN)
            MV__1 = [UN(jj),UT(jj)];
            FF(jj) = evalFcst_nlmpc(nlobj,X__, [MV__1;MV__], u_last, TRAJ__,[],nloptions);
        end
        figure(10)
        hold off
        contour(UN,UT,FF)
        hold on
        plot(u(1,i),u(2,i),'o')
        disp('DONE')
        pause
    end

    disp(info.ExitFlag)
    u_last = u(:,i);
    MVopt = info.MVopt;

    cost(i) = info.Cost;
    %%%%%%%%%%%%%%%% PLANT SIM

    x_dot_ = eval_model_continuous(x(:,i),u(:,i),slider);

    % Euler integration
    x(:,i+1) = x(:,i) + Ts*x_dot_;
    %                 toc;
    disp("NEXT ITERATION")
%     pause

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end
%             toc;

x_s = x(1,1:end-1); y_s = x(2,1:end-1); theta_s = x(3,1:end-1); S_p_x = -0.034*ones(size(theta_s)); S_p_y = x(4,1:end-1);
u_n = u(1,:);
u_t = u(2,:);

%% PLOT
set(0,'DefaultLineLineWidth',1.5);
figure,
time = time_sim_vec;
ax1 = subplot(3,2,1); plot(time,x_s), hold on, plot(time,traj(1,:)), xlabel('t [s]'), ylabel('x_S'),legend('x_S','xref_S'), subtitle('x_S tracking'), grid on
ax2 = subplot(3,2,3); plot(time,y_s), hold on, plot(time,traj(2,:)), xlabel('t [s]'), ylabel('y_S'),legend('y_S','yref_S'), subtitle('y_S tracking'), grid on
ax3 = subplot(3,2,5); plot(time,theta_s), hold on, plot(time,traj(3,:)),xlabel('t [s]'), ylabel('\theta_S'), legend('\theta_S','\thetaref_S'), subtitle('\theta_S tracking'), grid on
ax5 = subplot(3,2,2); plot(time,S_p_y),hold on, plot(time,traj(4,:)),xlabel('t [s]'), ylabel('S_ p_y'), legend('S_ p_y','ref'), subtitle('S_ p_y tracking'), grid on

figure
ax6 = subplot(2,1,1); plot(time,u_n), xlabel('t [s]'), ylabel('u_n'), subtitle("normal control"), grid on
ax7 = subplot(2,1,2); plot(time,u_t), xlabel('t [s]'), ylabel('u_t'), subtitle("tangential control"), grid on


linkaxes([ax1,ax2,ax3,ax5,ax6,ax7],'x');
xlim([ax1,ax2,ax3,ax5,ax6,ax7],[0 time(end)])


% helper.my_plot(params.t, traj, params.x_S, params.y_S, params.theta_S, params.S_p_y, params.u_n, params.u_t, controller.cost_function_vect, helper.convert_str2num(params.mode_vect));





