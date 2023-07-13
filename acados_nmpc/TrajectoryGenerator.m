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
        sample_time;
        set_plot = false;
    end

    methods
        function self = TrajectoryGenerator(T)
            % Constructor of the TrajectoryGenerator class
            % Input: T sample time [s]
            % Output: object
            self.sample_time = T;
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


        function [time, traj] = straight_line(self)
            time = self.t0:self.sample_time:self.tf;
            traj = ones(length(self.x0), length(time));
            for k = 1:length(time)
                s = quintic_(self,time(k))*norm(self.xf-self.x0);
                traj(:,k) = self.x0 + s*(self.xf-self.x0)/norm(self.xf-self.x0);
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
    end
end