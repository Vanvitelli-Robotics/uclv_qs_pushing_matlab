classdef helper
    properties(Constant)
        g = 9.81;
    end
    methods(Static)

        function R = my_rotz(theta)
            % Rotation matrix about z-axis (theta RAD)
            R = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];
        end

        function R = my_rotz_2d(theta)
            R_ = helper.my_rotz(theta);
            R = R_(1:2,1:2);
        end

        function mode_vect_num = convert_str2num(mode_vect_str)
            mode_vect_str = replace(mode_vect_str,"ST",'1');
            mode_vect_str = replace(mode_vect_str,"SL",'2');
            mode_vect_str = replace(mode_vect_str,"SR",'3');

            mode_vect_num = str2num(char(mode_vect_str));
        end

        function my_plot_robot(time, traj, x_s, y_s, theta_s, S_p_y, u_n, u_t)
            if length(traj(1,:)) < length(time)
                traj = [traj traj(:,end).*ones(length(traj(:,end)),length(time)-length(traj(1,:)))];
            end
            set(0,'DefaultLineLineWidth',1.5);
            figure,
            ax1 = subplot(3,2,1); plot(time,x_s), hold on, plot(time,traj(1,:)), xlabel('t [s]'), ylabel('x_S'),legend('x_S','xref_S'), subtitle('x_S tracking'), grid on
            ax2 = subplot(3,2,3); plot(time,y_s), hold on, plot(time,traj(2,:)), xlabel('t [s]'), ylabel('y_S'),legend('y_S','yref_S'), subtitle('y_S tracking'), grid on
            ax3 = subplot(3,2,5); plot(time,theta_s), hold on, plot(time,traj(3,:)),xlabel('t [s]'), ylabel('\theta_S'), legend('\theta_S','\thetaref_S'), subtitle('\theta_S tracking'), grid on
            %             ax4 = subplot(3,2,2); plot(time,S_p_x),hold on, plot(time,traj(4,:)),xlabel('t [s]'), ylabel('S_ p_x'), legend('S_ p_x','ref'), subtitle('S_ p_x tracking'), grid on
            ax5 = subplot(3,2,2); plot(time,S_p_y),hold on, plot(time,traj(4,:)),xlabel('t [s]'), ylabel('S_ p_y'), legend('S_ p_y','ref'), subtitle('S_ p_y tracking'), grid on

            figure
            ax6 = subplot(2,1,1); plot(time,u_n), xlabel('t [s]'), ylabel('u_n'), subtitle("normal control"), grid on
            ax7 = subplot(2,1,2); plot(time,u_t), xlabel('t [s]'), ylabel('u_t'), subtitle("tangential control"), grid on

            linkaxes([ax1,ax2,ax3,ax5,ax6,ax7],'x');
            xlim([ax1,ax2,ax3,ax5,ax6,ax7],[0 time(end)])
        end

        function my_plot(time, traj, x_s, y_s, theta_s, S_p_y, u_n, u_t,cost_function_vect, mode_vect)
            if length(traj(1,:)) < length(time)
                traj = [traj traj(:,end).*ones(length(traj(:,end)),length(time)-length(traj(1,:)))];
            end
            set(0,'DefaultLineLineWidth',1.5);
            figure,
            ax1 = subplot(3,2,1); plot(time,x_s), hold on, plot(time,traj(1,:)), xlabel('t [s]'), ylabel('x_S'),legend('x_S','xref_S'), subtitle('x_S tracking'), grid on
            ax2 = subplot(3,2,3); plot(time,y_s), hold on, plot(time,traj(2,:)), xlabel('t [s]'), ylabel('y_S'),legend('y_S','yref_S'), subtitle('y_S tracking'), grid on
            ax3 = subplot(3,2,5); plot(time,theta_s), hold on, plot(time,traj(3,:)),xlabel('t [s]'), ylabel('\theta_S'), legend('\theta_S','\thetaref_S'), subtitle('\theta_S tracking'), grid on
            %             ax4 = subplot(3,2,2); plot(time,S_p_x),hold on, plot(time,traj(4,:)),xlabel('t [s]'), ylabel('S_ p_x'), legend('S_ p_x','ref'), subtitle('S_ p_x tracking'), grid on
            ax5 = subplot(3,2,2); plot(time,S_p_y),hold on, plot(time,traj(4,:)),xlabel('t [s]'), ylabel('S_ p_y'), legend('S_ p_y','ref'), subtitle('S_ p_y tracking'), grid on

            figure
            ax6 = subplot(2,1,1); plot(time,u_n), xlabel('t [s]'), ylabel('u_n'), subtitle("normal control"), grid on
            ax7 = subplot(2,1,2); plot(time,u_t), xlabel('t [s]'), ylabel('u_t'), subtitle("tangential control"), grid on

            figure
            ax8 = subplot(2,1,1); plot(time, cost_function_vect), xlabel('t [s]'), ylabel('cost_function'), subtitle("Cost function"), grid on
            ax9 = subplot(2,1,2); plot(time(mode_vect==1), mode_vect(mode_vect==1),'*', "Color",'r'), hold on
            plot(time(mode_vect==2), mode_vect(mode_vect==2),'*', "Color",'g'),
            plot(time(mode_vect==3), mode_vect(mode_vect==3),'*', "Color",'y'),
            hold off
            xlabel('t [s]'), ylabel('mode'), subtitle("Mode"), legend('sticking','sliding left', 'sliding right'), grid on

            linkaxes([ax1,ax2,ax3,ax5,ax6,ax7,ax8,ax9],'x');
            xlim([ax1,ax2,ax3,ax5,ax6,ax7,ax8,ax9],[0 time(end)])
        end

        function my_animate(x, y, theta, rx, ry, t, step_time, traj)
            % Function to animate the trajectory tracking of the pusher slider system
            % Input: 2D position [x,y,theta]

            % Slider frame
            %     theta = -theta;
            p0_s = [x(1) y(1) 0];
            quat0_s = quaternion(helper.my_rotz(theta(1)),'rotmat','frame');
            figure
            slider_p = poseplot(quat0_s,p0_s, "MeshFileName","cad_models/cad_santal_centered_scaled_rotated_reduced.stl", "PatchFaceAlpha", 0.2, "ScaleFactor",0.001);%,PatchFaceColor="yellow");
            hold on

            p0_ref = [traj(1:2,1)' 0];
            quat0_ref = quaternion(helper.my_rotz(traj(3,1)),'rotmat','frame');
            slider_ref = poseplot(quat0_ref,p0_ref, "MeshFileName","cad_models/cad_santal_centered_scaled_rotated_reduced.stl", "PatchFaceAlpha", 0.2,"ScaleFactor",0.001);%,PatchFaceColor="yellow");


            % Pusher frame
            %     p0_p = helper.my_rotz(theta(1))*[rx(1) ry(1) 0]';
            %     quat0_p = quaternion(helper.my_rotz(0),'rotmat','frame');
            %     pusher_p = poseplot(quat0_p,p0_p, MeshFileName="cad model pusher\cad_model_pusher_1.stl",ScaleFactor=0.001);
            radius_pusher = 0.001;
            pusher_position = [x(1) y(1) 0]' + helper.my_rotz(theta(1))*[rx(1) ry(1) 0]';
            h = viscircles([pusher_position(1) pusher_position(2)],radius_pusher,"Color",'black',"LineWidth",6);
            viscircles([x(1) y(1)],radius_pusher,"Color",'blue',"LineWidth",6);

            hold off

            ax1 = gca; ax1.View = [-0.1303  -90];
            %             xlim(ax1,[-0.1 0.5])
            %             ylim(ax1,[-0.3 0.3])
            xlabel("x [m]")
            ylabel("y [m]")
            zlabel("z [m]")
            for t_k = 0:step_time:t(end)
                index_ = find(t>=t_k);
                i = index_(1);
                %             end
                %             for i = 2:1:length(x)
                h.Visible = false;
                position_s = [x(i) y(i) 0]';
                quat_s = quaternion([0 0 theta(i)],'rotvec');
                set(slider_p,"Orientation",quat_s,"Position",position_s);

                position_ref = [traj(1:2,i)' 0]';
                quat_ref = quaternion([0 0 traj(3,i)],'rotvec');
                set(slider_ref,"Orientation",quat_ref,"Position",position_ref);

                hold on
                pusher_position = [x(i) y(i) 0]'+ helper.my_rotz(theta(i))*[rx(i) ry(i) 0]';
                viscircles(ax1,[x(i) y(i)],radius_pusher,"Color",'blue',"LineWidth",6);
                h = viscircles(ax1,[pusher_position(1) pusher_position(2)],radius_pusher,"Color",'black',"LineWidth",6);
                plot(ax1,traj(1,:), traj(2,:), '-.','LineWidth',3,'Color','red');
                %                 pause(sampling_time)
                drawnow
            end
        end

        function [x_s, y_s, theta_s, S_p_x, S_p_y, u_n, u_t, time_sim_vec] = open_loop_matlab(plant,x0,u_n,u_t, time_sim,sample_time, sim_noise)
            % OPEN LOOP SIMULATION for variable shape

            % Time of overall simulation
            time_sim_vec = 0:sample_time:time_sim;
            time_sim_ = length(time_sim_vec);

            % State
            x = zeros(plant.nx, time_sim_+1);
            x(:,1) = x0;

            % Input Signal
            u = repmat([u_n;u_t],1, time_sim_);

            % Delay buffer to simulate the delay of the plant
            delay_buff_plant = ceil(plant.time_delay / sample_time);
            u_buff_plant = zeros(plant.nu,delay_buff_plant);

            tic;
            for i = 1:time_sim_

                % noise simulation
                if(sim_noise == true)
                    x(:,i) = x(:,i) + [1e-5*randn(1,2) randn()*1e-3 randn()*1e-4 randn()*1e-4]';
                end

                if delay_buff_plant == 0
                    x_dot_ = plant.evalModelVariableShape(x(:,i),u(:,i)); %plant.eval_model(x(:,i),u(:,i));
%                     x_dot_ = plant.eval_model(x(:,i),u(:,i));
%                     x_dot_ = plant.eval_model_variable_shape(x(:,i),u(:,i));
                else
                    x_dot_ = plant.evalModelVariableShape(x(:,i),u_buff_plant(:,end));
%                     x_dot_ = plant.eval_model(x(:,i),u_buff_plant(:,end));
                    u_buff_plant = [u(:,i) u_buff_plant(:,1:end-1)];
                end

                % Euler integration
                x(:,i+1) = x(:,i) + sample_time*x_dot_;

            end
            toc;


            x_s = x(1,1:end-1); y_s = x(2,1:end-1); theta_s = x(3,1:end-1); %S_p_x = x(4,1:end-1); S_p_y = x(5,1:end-1);
            S_p = plant.SP.evalSpline(plant.SP.FC,x(4,1:end-1));
            S_p_x = S_p(:,1)';
            S_p_y = S_p(:,2)';
            u_n = u(1,:);
            u_t = u(2,:);
        end

        function [x_s, y_s, theta_s, S_p_x, S_p_y, u_n, u_t, time_sim_vec,mode_vect, found_sol] = closed_loop_matlab(plant, controller,x0, time_sim, print_, sim_noise, debug_cost,disturbance_)
            % CLOSED LOOP SIMULATION
            index = 0;
            % Time of overall simulation
            time_sim_vec = 0:controller.sample_time:time_sim;
            time_sim_ = length(time_sim_vec);

            % State and control variables
            x = zeros(plant.nx, time_sim_+1);
            x(:,1) = x0;
            u = zeros(plant.nu, time_sim_);
            mode_vect = string(zeros(time_sim_,1));
            found_sol = false(length(time_sim_vec),1);

            % Delay buffer to simulate the delay of the plant
            delay_buff_plant = ceil(plant.time_delay / controller.sample_time);
            u_buff_plant = zeros(plant.nu,delay_buff_plant);


            %             tic;
            for i = 1:time_sim_
                disp("Time: " + i*controller.sample_time);
                if disturbance_ == true
                    if i == 33
                        disp("Disturbance")
                        x(2,i) = x(2,i) + 0.02;
                        x(4,i) = x(4,i) - 0.02;
                    end
                end


                % noise simulation
                if(sim_noise == true)
                    x(:,i) = x(:,i) + [1e-5*randn(1,2) randn()*1e-3 randn()*1e-4]';
                end

                xk_sim = controller.delay_buffer_sim(plant, x(:,i));

                % solve OCP
                tic;
                u(:,i) = controller.solve(xk_sim,i+controller.delay_buff_comp);
                toc;

                controller.u_buff_contr = [u(:,i) controller.u_buff_contr(:,1:end-1)];
                status = controller.ocp_solver.get('status');
                if status~=0
                    disp('acados ocp solver failed');
                    found_sol(i) = false;
                    %                         keyboard
                else
                    found_sol(i) = true;
                end
            

                if print_ == true
                    controller.ocp_solver.print;
                    
                    sqp_iter = controller.ocp_solver.get('sqp_iter');
                    time_tot = controller.ocp_solver.get('time_tot');
                    time_lin = controller.ocp_solver.get('time_lin');
                    time_qp_sol = controller.ocp_solver.get('time_qp_sol');

                    fprintf('\nstatus = %d, sqp_iter = %d, time_int = %f [ms] (time_lin = %f [ms], time_qp_sol = %f [ms])\n',...
                        status, sqp_iter, time_tot*1e3, time_lin*1e3, time_qp_sol*1e3);
                end

                

                if debug_cost == true
                    if i*controller.sample_time < 0.06 || (i*controller.sample_time > 0.85 && i*controller.sample_time < 1.2)% && mod(i,5)==0)
                        index = index + 1;
                        u_min = helper.debug_cost_function(x(:,i), controller.u_n_lb, controller.u_t_lb, controller.u_n_ub, controller.u_t_ub, controller, plant, index, i*controller.sample_time);
                        %                         u(:,i) = u_min;
                    end
                    if index == 9
                        index = 0;
                        figure;
                    end
                end


                %%%%%%%%%%%%%%%% PLANT SIM

                if delay_buff_plant == 0
%                     [x_dot_, mode_] = plant.eval_model(x(:,i),u(:,i));
                    [x_dot_, mode_] = plant.eval_model_variable_shape(x(:,i),u(:,i));
                else
%                     [x_dot_, mode_] = plant.eval_model(x(:,i),u_buff_plant(:,end));
                    [x_dot_, mode_] = plant.eval_model_variable_shape(x(:,i),u_buff_plant(:,end));

                    u_buff_plant = [u(:,i) u_buff_plant(:,1:end-1)];
                end

                mode_vect(i) = mode_;



                % Euler integration
                x(:,i+1) = x(:,i) + controller.sample_time*x_dot_;
                %                 toc;

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            end
            %             toc;

            x_s = x(1,1:end-1); y_s = x(2,1:end-1); theta_s = x(3,1:end-1); %S_p_x = -0.034*ones(size(theta_s)); S_p_y = x(4,1:end-1);
            S_p = plant.SP.evalSpline(plant.SP.FC,x(4,1:end-1));
            S_p_x = S_p(:,1)';
            S_p_y = S_p(:,2)';
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

        function params = save_parameters(name_exp, x, u, t, mode_vect, params)
            %             params = struct;

            params.t = t;
            params.x_S = x(1,:);
            params.y_S = x(2,:);
            params.theta_S = x(3,:);
            %                         params.S_p_x = x(4,:);
            params.S_p_y = x(4,:);
            params.u_n = u(1,:);
            params.u_t = u(2,:);
            if nargin > 4
                params.mode_vect = mode_vect;
            end
            name_exp_ext = strcat(name_exp,'.mat');
            save(name_exp_ext,'params');
        end

        function cost_fcn_ = evaluate_cost_function(xn, controller, plant,t, u_vect)
            cost_fcn_ = 0;
            time_index = round(t/controller.sample_time);
            for n = 1:controller.Hp
                y_ref = controller.get_y_ref(time_index+n);
                cost_fcn_ = cost_fcn_ + (xn-y_ref(1:4))'*controller.W_x*(xn-y_ref(1:4))+(u_vect(:,n)-y_ref(5:end))'*controller.W_u*(u_vect(:,n)-y_ref(5:end));            %f_cost_k(-xn_ref+xn,u_vect(:,n));
                x_dot = plant.eval_model(xn,u_vect(:,n));
                xn = xn + controller.sample_time*x_dot;
            end
            yn_ref = controller.get_y_ref(time_index+controller.Hp);
            cost_fcn_ = cost_fcn_ + (xn-yn_ref(1:4))'*controller.W_x_e*(xn-yn_ref(1:4));
        end

        function u_min = debug_cost_function(x0, u_n_lb, u_t_lb, u_n_ub, u_t_ub, controller, plant, index, t)
            import casadi.*
            %             f_cost_k = Function('f_k',{controller.sym_model.sym_x, controller.sym_model.sym_u},{controller.ocp_solver.model_struct.cost_expr_ext_cost});
            %             f_cost_ke = Function('f_ke',{controller.sym_model.sym_x},{controller.ocp_solver.model_struct.cost_expr_ext_cost_e});
            %     controller.clear_variables();
            %     u_opt = controller.solve(x0);
            %     u_vect = controller.ocp_solver.get('u');
            u_vect = controller.ocp_solver.get('u');
            u_opt = u_vect(:,1);
            %             cost_values = [];
            cost_fcn_ = 0;
            UN = u_n_lb:0.005:u_n_ub;
            UT= u_t_lb:0.005:u_t_ub;
            [X,Y] = meshgrid(UN,UT);


            cost_values = zeros(size(X));

            %
            %             if true || index <10
            %                 for i = 1:length(UN)
            %                     for j = 1:length(UT)
            %                         u_vect(:,1) = [UN(i);UT(j)];
            %                         xn = x0;
            %                         % Evaluate cost function
            %                         for n = 1:controller.Hp-1
            %                             xn_ref = controller.y_ref(1:5,controller.index_ref+n-1);
            %                             cost_fcn_ = cost_fcn_ + f_cost_k(-xn_ref+xn,u_vect(:,n));
            %                             x_dot = plant.eval_model(xn,u_vect(:,n));
            %                             xn = xn + controller.sample_time*x_dot;
            %                         end
            %                         xn_ref = controller.y_ref(1:5,controller.index_ref+controller.Hp-1);
            %                         cost_fcn_ = cost_fcn_ + f_cost_ke(-xn_ref+xn);
            %                         % %             cost_values = [cost_values; full(cost_fcn_)];
            %                         cost_values(j,i) = full(cost_fcn_);
            %                         cost_fcn_ = 0;
            %                     end
            %                 end
            %             end

            time_index = round(t/controller.sample_time);%,t,controller.Hp,(t/controller.sample_time)-round((t/controller.sample_time))

            if true || index <10
                for i = 1:numel(X)
                    u_vect(:,1) = [X(i);Y(i)];
                    xn = x0;
                    % Evaluate cost function
                    for n = 1:controller.Hp
                        y_ref = controller.get_y_ref(time_index+n);
                        cost_fcn_ = cost_fcn_ + 0.5*((xn-y_ref(1:4))'*controller.W_x*(xn-y_ref(1:4))+(u_vect(:,n)-y_ref(5:end))'*controller.W_u*(u_vect(:,n)-y_ref(5:end)));            %f_cost_k(-xn_ref+xn,u_vect(:,n));
                        x_dot = plant.eval_model(xn,u_vect(:,n));
                        xn = xn + controller.sample_time*x_dot;
                    end
                    yn_ref = controller.get_y_ref(time_index+controller.Hp);
                    cost_fcn_ = cost_fcn_ + 0.5*((xn-yn_ref(1:4))'*controller.W_x_e*(xn-yn_ref(1:4)));
                    % %             cost_values = [cost_values; full(cost_fcn_)];
                    cost_values(i) = (cost_fcn_);
                    cost_fcn_ = 0;
                end

            end
            [min_, ind_raw] = min(cost_values(:));
            u_min = [X(ind_raw); Y(ind_raw)];

            if true
                subplot(3,3,index), contour(X,Y,cost_values)
                %                 [min_, ind_raw] = min(cost_values);
                %                 [min_2, ind_col] = min(min_);


                xlabel('u_n [m/s]')
                ylabel('u_t [m/s]')
                subtitle(strcat("Time [s]: ",num2str(t)))
                hold on
                plot(u_opt(1), u_opt(2),'*')
                %                             plot(u_vect(1,2), u_vect(2,2),'*')
                %                             plot(u_vect(1,3), u_vect(2,3),'*')
                hold off
                pause
            end


        end
    end


end