classdef TrajectoryGenerator < handle
    % This class creates the symbolic nonlinear pusher-slider model by using acados.
    % It can be used also to simulate the plant.
    %
    % PusherSliderModel Properties:
    %    slider_params -
    %
    % PusherSliderModel Methods:
    %    PusherSliderModel - Constructor

    properties
        x0; % Start point of the desired trajectory [m]
        xf; % Final point of the desired trajectory [m]
        tf; % Time to execute the overall trajectory [s]
        t0;
        waypoints_;
        waypoints_velocities;
        sample_time;
        set_plot = false;
        vel;
    end

    methods
        function self = TrajectoryGenerator(T, vel)
            % Constructor of the TrajectoryGenerator class
            % Input: T sample time [s]
            % Output: object
            self.sample_time = T;
            self.vel = vel;
        end

        function set_target(self,x0,xf,t0,tf)
            self.x0 = x0;
            self.xf = xf;
            self.tf = tf;
            self.t0 = t0;
        end

        function s = quintic_(self, time)
            tau = time/self.tf;
            s = 6*tau.^5-15*tau.^4+10*tau.^3;
        end

        function [time, traj] = straight_line(self, auto_angle)
            time = self.t0:self.sample_time:self.tf;
            traj = ones(length(self.x0), length(time));

            for k = 1:length(time)
                s = quintic_(self,time(k))*norm(self.xf-self.x0);
                traj(:,k) = self.x0 + s*(self.xf-self.x0)/norm(self.xf-self.x0);
            end

            % Trajectory for angle -> it has to be faster than the
            % trajectory of x and y
            if auto_angle == true
                tf_angle = self.tf/2;
                traj_angle = ones(1, length(time));
                time_angle = self.t0:self.sample_time:tf_angle;
                for j = 1:length(time_angle)
                    s = quintic_(self,time_angle(j))*norm(self.xf(3)-self.x0(3));
                    traj_angle(j) = self.x0(3) + s*(self.xf(3)-self.x0(3))/norm(self.xf(3)-self.x0(3));
                end
                traj_angle(length(time_angle):length(time)) = traj_angle(3,end).*ones(1,length(time)-length(time_angle));
                traj = [traj(1:2,:); traj_angle; traj(4:5,:)];
            end

            if self.set_plot == true
                figure, plot(time,traj), grid on
                legend
                xlabel("time [s]")
                ylabel("traj [m]")
                if length(self.xf) > 1
                    figure, plot(traj(1,:), traj(2,:),'*'), grid on
                    xlabel("x [m]")
                    ylabel("y [m]")
                end
            end

        end

        function [time, traj] = waypoint_gen_fixed_angle(self)
            delta_p = (self.waypoints_(2:end,:)-self.waypoints_(1:end-1,:));
            times = vecnorm(delta_p')./self.waypoints_velocities(:)';

            time = 0;
            traj = [self.waypoints_(1,1:2)';self.x0(3:4); 0];
            for i=1:1:length(self.waypoints_)-1
                time_i = (time(end)+self.sample_time):self.sample_time:(time(end)+times(i));
                time = [time time_i];
                x_i = linspace(self.waypoints_(i,1),self.waypoints_(i+1,1),length(time_i));
                y_i = linspace(self.waypoints_(i,2),self.waypoints_(i+1,2),length(time_i));
                traj = [traj, [x_i;y_i;repmat([self.x0(3:4);0],1,length(x_i))]];   
            end 

        end 
        function [time, traj] = waypoints_gen(self)
            Fs = 1/self.sample_time;

            delta_p = abs(self.waypoints_(2:end,:)-self.waypoints_(1:end-1,:));
            
%             times = vecnorm(delta_p')/self.vel;
            times = vecnorm(delta_p')./self.waypoints_velocities(:)';
            times = [0 times];
            time_tf = cumsum(times);

            % create trajectory
%             eulerAngs = zeros(size(self.waypoints_,1),3);
%             eulerAngs(1,:) = [self.x0(3) 0 0];
%             eulerAngs(2:end,1) = atan2(self.waypoints_(2:end,2)-self.waypoints_(1:end-1,2),self.waypoints_(2:end,1)-self.waypoints_(1:end-1,1));
%             q = quaternion(eulerAngs,"euler","ZYX","frame");

            orientations = zeros(3,3,size(self.waypoints_,1));
            orientations(:,:,1) = eye(3);
            for i=2:size(self.waypoints_,1)-1
                ang = atan2(self.waypoints_(i+1,2)-self.waypoints_(i,2),self.waypoints_(i+1,1)-self.waypoints_(i,1));
                orientations(:,:,i) = helper.my_rotz(-ang);
            end
            orientations(:,:,end) = eye(3);


            trajectory = waypointTrajectory(self.waypoints_, time_tf, 'SampleRate', Fs);% 'AutoPitch', true);% 'Orientation',orientations);

            % lookup pose information for entire trajectory
            [pos, orient] = lookupPose(trajectory, time_tf(1):1/Fs:time_tf(end));

            % Plot generated positions and specified waypoints.
            if self.set_plot == true
                plot(pos(:,1),pos(:,2), self.waypoints_(:,1),self.waypoints_(:,2), '--o')
                title('Position')
                xlabel('X (m)')
                ylabel('Y (m)')
                zlabel('Z (m)')
                legend({'Position', 'Waypoints'})
            end

%             yaw = quat2angle(orient);
            eulerAngles = quat2eul(orient);
            yaw = eulerAngles(:,1);

            traj = [pos(:,1)'; pos(:,2)'; yaw'; self.x0(4)*ones(1,length(pos(:,1))); zeros(1,length(pos(:,1)))];
            tInfo = waypointInfo(trajectory);
            time = 0:(1/trajectory.SampleRate):tInfo.TimeOfArrival(end);
        end
    end
end