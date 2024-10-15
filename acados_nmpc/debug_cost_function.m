
function debug_cost_function(x0, u_n_lb, u_t_lb, u_n_ub, u_t_ub, controller, plant, index, t)
    import casadi.*
    f_cost_k = Function('f_k',{controller.sym_model.sym_x, controller.sym_model.sym_u},{controller.ocp_model.model_struct.cost_expr_ext_cost});
    f_cost_ke = Function('f_ke',{controller.sym_model.sym_x},{controller.ocp_model.model_struct.cost_expr_ext_cost_e});
    %     controller.clear_variables();
    %     u_opt = controller.solve(x0);
    %     u_vect = controller.ocp_solver.get('u');
    u_vect = controller.ocp_solver.get('u');
    u_opt = u_vect(:,1);
    cost_values = [];
    cost_fcn_ = 0;
    UN = u_n_lb:0.005:u_n_ub;
    UT= u_t_lb:0.005:u_t_ub;
    [X,Y] = meshgrid(UN,UT);


    if index <10
        for i = 1:length(UN)
            for j = 1:length(UT)
                u_vect(:,1) = [UN(i);UT(j)];
                xn = x0;
                % Evaluate cost function
                for n = 1:controller.Hp-1
                    xn_ref = controller.y_ref(1:5,controller.index_ref+n-1);
                    cost_fcn_ = cost_fcn_ + f_cost_k(-xn_ref+xn,u_vect(:,n));
                    x_dot = plant.eval_model(xn,u_vect(:,n));
                    xn = xn + controller.sample_time*x_dot;
                end
                xn_ref = controller.y_ref(1:5,controller.index_ref+controller.Hp-1);
                cost_fcn_ = cost_fcn_ + f_cost_ke(-xn_ref+xn);
                % %             cost_values = [cost_values; full(cost_fcn_)];
                cost_values(j,i) = full(cost_fcn_);
                cost_fcn_ = 0;
            end
        end


        if index == 11
            figure
        end
        subplot(3,3,index), contour(X,Y,cost_values)
        xlabel('u_n [m/s]')
        ylabel('u_t [m/s]')
        subtitle(strcat("Time [s]: ",num2str(t)))
        hold on
        plot(u_opt(1), u_opt(2),'*')
        hold off
    end

end