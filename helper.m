classdef helper
    properties(Constant)
        g = 9.81;
    end
    methods(Static)
        function integral = DoubleGaussQuad(fun1,a,b,c,d)
            %Change of Variables
            h1 = (b-a)/2;
            h2 = (b+a)/2;
            h3 = (d-c)/2;
            h4 = (d+c)/2;

            %Define Weights (for 3 points)
            w1 = 1;
            w2 = 1;

            %Define Points
            x1 = sqrt(1/3);
            x2 = -sqrt(1/3);

            integral = h1 * h3 * (w1*w1*fun1(h1*x1+h2, h1*x1+h2) + w1*w2*fun1(h1*x1+h2, h1*x2+h2) +...
                w2*w1*fun1(h1*x2+h2, h1*x1+h2) + w2*w2*fun1(h1*x2+h2, h1*x2+h2) );
        end

        function tau_max = tau_max_func(mu_sg, m, g, area, xwidth, ywidth)
            n_f_integrand = @(p1, p2) (mu_sg * m * g / area) * sqrt([p1; p2; 0]' * [p1; p2; 0]);
            tau_max = helper.DoubleGaussQuad(n_f_integrand, -xwidth / 2, xwidth / 2, -ywidth / 2, ywidth / 2);
        end

        function R = my_rotz(theta)
            % Rotation matrix about z-axis (theta RAD)
            R = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];
        end

        function my_plot(time, traj, x_s, y_s, theta_s, S_p_x, S_p_y, u_n, u_t)
            if length(traj(1,:)) < length(time)
                traj = [traj traj(:,end).*ones(length(traj(:,end)),length(time)-length(traj(1,:)))];
            end
            set(0,'DefaultLineLineWidth',1.5);
            figure
            subplot(3,2,1), plot(time,x_s), hold on, plot(time,traj(1,:)), xlabel('t [s]'), ylabel('x_S'),legend('x_S','xref_S'), subtitle('x_S tracking'), grid on
            subplot(3,2,2), plot(time,y_s), hold on, plot(time,traj(2,:)), xlabel('t [s]'), ylabel('y_S'),legend('y_S','yref_S'), subtitle('y_S tracking'), grid on
            subplot(3,2,3), plot(time,theta_s), hold on, plot(time,traj(3,:)),xlabel('t [s]'), ylabel('\theta_S'), legend('\theta_S','\thetaref_S'), subtitle('\theta_S tracking'), grid on
            subplot(3,2,4), plot(time,S_p_x),hold on, plot(time,traj(4,:)),xlabel('t [s]'), ylabel('S_ p_x'), legend('S_ p_x','ref'), subtitle('S_ p_x tracking'), grid on
            subplot(3,2,5), plot(time,S_p_y),hold on, plot(time,traj(5,:)),xlabel('t [s]'), ylabel('S_ p_y'), legend('S_ p_y','ref'), subtitle('S_ p_y tracking'), grid on

            figure
            subplot(2,1,1), plot(time,u_n), xlabel('t [s]'), ylabel('u_n'), subtitle("normal control"), grid on
            subplot(2,1,2), plot(time,u_t), xlabel('t [s]'), ylabel('u_t'), subtitle("tangential control"), grid on
        end


        function my_animate(x, y, theta, rx, ry,sampling_time, traj)
            % Function to animate the trajectory tracking of the pusher slider system
            % Input: 2D position [x,y,theta]

            % Slider frame
            %     theta = -theta;
            p0_s = [x(1) y(1) 0];
            quat0_s = quaternion(helper.my_rotz(theta(1)),'rotmat','frame');
            slider_p = poseplot(quat0_s,p0_s, MeshFileName="cad_models/cuboide_santal.stl", ScaleFactor=0.001);%,PatchFaceColor="yellow");
            hold on

            % Pusher frame
            %     p0_p = helper.my_rotz(theta(1))*[rx(1) ry(1) 0]';
            %     quat0_p = quaternion(helper.my_rotz(0),'rotmat','frame');
            %     pusher_p = poseplot(quat0_p,p0_p, MeshFileName="cad model pusher\cad_model_pusher_1.stl",ScaleFactor=0.001);
            radius_pusher = 0.001;
            pusher_position = [x(1) y(1) 0]' + helper.my_rotz(theta(1))*[rx(1) ry(1) 0]';
            viscircles([pusher_position(1) pusher_position(2)],radius_pusher,"Color",'red',"LineWidth",6);

            hold off

            ax1 = gca; ax1.View = [-0.1303  -90];
            xlim(ax1,[-0.1 0.5])
            ylim(ax1,[-0.3 0.3])
            xlabel("x [m]")
            ylabel("y [m]")
            zlabel("z [m]")

            for i = 2:1:length(x)
                position_s = [x(i) y(i) 0]';
                quat_s = quaternion([0 0 theta(i)],'rotvec');
                set(slider_p,Orientation=quat_s,Position=position_s);
                hold on
                pusher_position = [x(i) y(i) 0]'+ helper.my_rotz(theta(i))*[rx(i) ry(i) 0]';
                viscircles([pusher_position(1) pusher_position(2)],radius_pusher,"Color",'red',"LineWidth",6);
                plot(traj(1,:), traj(2,:), '-.','LineWidth',3,'Color','b');
                pause(sampling_time)
                drawnow
            end
        end

    end

end