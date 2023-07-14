classdef helper
    properties(Constant)
        g = 9.81;
    end
    methods(Static)

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
            ax1 = subplot(3,2,1); plot(time,x_s), hold on, plot(time,traj(1,:)), xlabel('t [s]'), ylabel('x_S'),legend('x_S','xref_S'), subtitle('x_S tracking'), grid on
            ax2 = subplot(3,2,3); plot(time,y_s), hold on, plot(time,traj(2,:)), xlabel('t [s]'), ylabel('y_S'),legend('y_S','yref_S'), subtitle('y_S tracking'), grid on
            ax3 = subplot(3,2,5); plot(time,theta_s), hold on, plot(time,traj(3,:)),xlabel('t [s]'), ylabel('\theta_S'), legend('\theta_S','\thetaref_S'), subtitle('\theta_S tracking'), grid on
            ax4 = subplot(3,2,2); plot(time,S_p_x),hold on, plot(time,traj(4,:)),xlabel('t [s]'), ylabel('S_ p_x'), legend('S_ p_x','ref'), subtitle('S_ p_x tracking'), grid on
            ax5 = subplot(3,2,4); plot(time,S_p_y),hold on, plot(time,traj(5,:)),xlabel('t [s]'), ylabel('S_ p_y'), legend('S_ p_y','ref'), subtitle('S_ p_y tracking'), grid on

            figure
            ax6 = subplot(2,1,1); plot(time,u_n), xlabel('t [s]'), ylabel('u_n'), subtitle("normal control"), grid on
            ax7 = subplot(2,1,2); plot(time,u_t), xlabel('t [s]'), ylabel('u_t'), subtitle("tangential control"), grid on
            linkaxes([ax1,ax2,ax3,ax4,ax5,ax6,ax7],'x');
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

        function [x_s, y_s, theta_s, S_p_x, S_p_y, u_n, u_t, time_sim_vec] = closed_loop_matlab(plant, controller,x0, time_sim, print_, Wx1, Wx2,Wxe,Wu)
            % CLOSED LOOP SIMULATION

            % Time of overall simulation
            time_sim_vec = 0:controller.sample_time:time_sim;
            time_sim_ = length(time_sim_vec);

            % State and control variables
            x = zeros(plant.nx, time_sim_+1);
            x(:,1) = x0;
            u = zeros(plant.nu, time_sim_);

            % Delay buffer to simulate the delay of the plant
            delay_buff_plant = ceil(plant.time_delay / controller.sample_time);
            u_buff_plant = zeros(plant.nu,delay_buff_plant);

            % Delay buffer to compensate the delay with the controller
            delay_buff_comp = ceil(controller.delay_compensation/controller.sample_time);
            u_buff_contr = zeros(plant.nu, delay_buff_comp);

            tic;
            for i = 1:time_sim_
                if i < time_sim_/5
                    controller.update_cost_function(Wx1,Wu,Wxe,1,controller.Hp-1);
                else
                    controller.update_cost_function(Wx2,Wu,Wxe,1,controller.Hp-1);
                end

                xk_sim = x(:,i);
                for k = 1 : delay_buff_comp
                    x_dot_sim = plant.eval_model(xk_sim,u_buff_contr(:,end-k+1));
                    x_sim = xk_sim + controller.sample_time*x_dot_sim;
                    xk_sim = x_sim;
                end

                % solve OCP
                u(:,i) = controller.solve(xk_sim);
                u_buff_contr = [u(:,i) u_buff_contr(:,1:end-1)];

                if print_ == true
                    status = controller.ocp_solver.get('status');
                    sqp_iter = controller.ocp_solver.get('sqp_iter');
                    time_tot = controller.ocp_solver.get('time_tot');
                    time_lin = controller.ocp_solver.get('time_lin');
                    time_qp_sol = controller.ocp_solver.get('time_qp_sol');

                    fprintf('\nstatus = %d, sqp_iter = %d, time_int = %f [ms] (time_lin = %f [ms], time_qp_sol = %f [ms])\n',...
                        status, sqp_iter, time_tot*1e3, time_lin*1e3, time_qp_sol*1e3);
                    if status~=0
                        disp('acados ocp solver failed');
                        keyboard
                    end
                end

                if delay_buff_plant == 0
                    x_dot_ = plant.eval_model(x(:,i),u(:,i));
                else
                    x_dot_ = plant.eval_model(x(:,i),u_buff_plant(:,end));
                    u_buff_plant = [u(:,i) u_buff_plant(:,1:end-1)];
                end

                % Euler integration
                x(:,i+1) = x(:,i) + controller.sample_time*x_dot_;

            end
            toc;

            x_s = x(1,1:end-1); y_s = x(2,1:end-1); theta_s = x(3,1:end-1); S_p_x = x(4,1:end-1); S_p_y = x(5,1:end-1);
            u_n = u(1,:);
            u_t = u(2,:);
        end

        function [x_s, y_s, theta_s, S_p_x, S_p_y, u_n, u_t, time_sim_vec] = closed_loop_simulink(time_sim)
            %             hws = get_param(model, 'modelworkspace');
            sim('simulation_model_closed_loop',time_sim);
            x_s = out.signals.values(:,1);
            y_s = out.signals.values(:,2);
            theta_s = out.signals.values(:,3);
            S_p_x = out.signals.values(:,4);
            S_p_y = out.signals.values(:,5);
            u_n = u_control.signals.values(:,1);
            u_t = u_control.signals.values(:,2);
            time_sim_vec = out.time;
        end

    end


end