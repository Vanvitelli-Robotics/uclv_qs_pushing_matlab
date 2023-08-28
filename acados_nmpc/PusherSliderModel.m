classdef PusherSliderModel < casadi.Callback
    % This class creates the symbolic nonlinear pusher-slider model by using acados.
    % It can be used also to simulate the plant.
    %
    % PusherSliderModel Properties:
    %    slider_params - Struct of the slider params, containing:
    %                       % slider.mu_sg   friction coefficient between slider and ground
    %                       % slider.mu_sp   friction coefficient between slider and pusher
    %                       % slider.ywidth  width of the slider along x-direction [m]
    %                       % slider.xwidth  width of the slider along y-direction [m]
    %                       % slider.area    slider area [m^2]
    %                       % slider.m       slider mass [kg]
    %    nx - Number of the state variables
    %    nu - Number of the control variables
    %    sym_model - Struct to create the symbolic model
    %
    % PusherSliderModel Methods:
    %    PusherSliderModel - Constructor
    %    contact_surface - Check if there is contact between pusher and slider
    %    eval_model - Evaluate numerical nonlinear model for simulation
    %    symbolic_model - Evaluate symbolic nonlinear model

    properties
        slider_params = struct;

        nx = 4;     % number of state variables
        nu = 2;     % number of inputs

        sym_model = struct;
        time_delay;
        xdot_func; % symbolic expression of model

        cad_model = struct;
        has_cad_model = false

        SP; % Spline object interpolating cad model points

    end

    methods

        % Constructor
        function self = PusherSliderModel(name, slider_parameters, time_delay,cad_model_path,order_spline,z_limit)
            % Constructor of the pusher_slider model.
            % Input: name = string
            %        slider_parameters = struct
            % Ouptut: object instance of PusherSliderModel
            self@casadi.Callback();
            self.slider_params=slider_parameters;
            self.slider_params.f_max = self.slider_params.mu_sg*self.slider_params.m*helper.g;
            self.slider_params.tau_max = self.tau_max_func(self.slider_params.mu_sg, self.slider_params.m, helper.g, self.slider_params.area, self.slider_params.xwidth, self.slider_params.ywidth);
            self.slider_params.c_ellipse = self.slider_params.tau_max/self.slider_params.f_max;
            self.set_delay(time_delay);
            self.open_cad_model(cad_model_path,order_spline,z_limit);
            construct(self, name);
        end

        function open_cad_model(self, cad_model_path,p,z_limit)
            % This function open and save the cad model of the Slider
            if(not(isempty(cad_model_path)))
                self.cad_model.path = cad_model_path;
                self.cad_model.stl = stlread(self.cad_model.path);
                self.cad_model.DT = delaunayTriangulation(self.cad_model.stl.Points);
                [self.cad_model.T,self.cad_model.Xb] = freeBoundary(self.cad_model.DT);
                self.cad_model.TR = triangulation(self.cad_model.T,self.cad_model.Xb);
                self.cad_model.TriangleCenters = incenter(self.cad_model.TR);
                self.cad_model.TriangleNormals = faceNormal(self.cad_model.TR);
                self.cad_model.scale_factor = 1000; % mesh scale
                self.has_cad_model = true;
                disp("Cad model saved");

                self.getSpline(p,z_limit);

            else
                self.has_cad_model = false;
                disp("No cad model specified");
            end
        end

        function Pxy_sorted = sortCadPoints(self,z_limit)
            PC = pcread('../cad_models/santal_planar_surface_simplified.ply');
            Pxy = PC.Location(:,1:2);
%             P = self.cad_model.stl.Points;
%             Pxy = P(abs(P(:,3))<z_limit,:);
            Pxy = Pxy(:,1:2);

            [~, ind] = min(Pxy(:,1));
            Pxy_tmp = Pxy(ind,:);
            Pxy(ind,:) = Inf;

            Pxy_sorted = zeros(size(Pxy));
            Pxy_sorted(1,:) = Pxy_tmp;

            for i = 2:length(Pxy)
                [~, ind_tmp] = min(vecnorm((Pxy-Pxy_tmp)'));
                Pxy_tmp = Pxy(ind_tmp,:);
                Pxy_sorted(i,:) = Pxy_tmp;
                Pxy(ind_tmp,:) = Inf;
            end

            Pxy_sorted = Pxy_sorted.*(1/self.cad_model.scale_factor);
            Pxy_sorted(end+1,:) = Pxy_sorted(1,:);
            Pxy_sorted = flipud(Pxy_sorted);
        end

        function getSpline(self,p,z_limit)
            % Points of contact P
            P = self.sortCadPoints(z_limit);
            
            n = length(P);

            % Knots vector S = [s0 s1, ..., sm]
            m = n+p+1-2*p;
            a = 0; b = sum(vecnorm(diff(P)'));
            S_ = linspace(a,b,m);
            S = [a*ones(1,p) S_ b*ones(1,p)];

            self.SP = bspline_shape(S,P,p);
            self.SP.getSymbolicSpline(p);
            self.SP.getSymboliSplineDot(p);
            self.SP.getNormalTangentialVersors;
        end

        function print_cad_model(self)
            figure,
            trisurf(self.cad_model.T,self.cad_model.Xb(:,1),self.cad_model.Xb(:,2),self.cad_model.Xb(:,3), ...
                'FaceColor','cyan','FaceAlpha',0.8);
            axis equal
            %
            %             hold on
            %             quiver3(self.cad_model.TriangleCenters(idx,1),self.cad_model.TriangleCenters(idx,2),self.cad_model.TriangleCenters(idx,3), ...
            %             normal(1),normal(2),0,100,'color','r');
            %             hold on
            %             quiver3(self.cad_model.TriangleCenters(idx,1),self.cad_model.TriangleCenters(idx,2),self.cad_model.TriangleCenters(idx,3), ...
            %                                 tangential(1),tangential(2),0,100,'color','b');
            %             plot(p(1)*1000,p(2)*1000,'ro')
        end

        function set_delay(self, time_delay)
            self.time_delay = time_delay;
        end

        function integral = DoubleGaussQuad(self,fun1,a,b,c,d)
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

        function tau_max = tau_max_func(self, mu_sg, m, g, area, xwidth, ywidth)
            n_f_integrand = @(p1, p2) (mu_sg * m * g / area) * sqrt([p1; p2; 0]' * [p1; p2; 0]);
            tau_max = self.DoubleGaussQuad(n_f_integrand, -xwidth / 2, xwidth / 2, -ywidth / 2, ywidth / 2);
        end

        % Contact surface
        function contact = contact_surface(self,S_p_x,S_p_y)
            % This method evaluates if the pusher at (S_p_x,S_p_y) is in contact with the slider.
            % Input: S_p_x  x-coordinate of the pusher position w.r.t slider frame
            %        S_p_y  y-coordinate of the pusher position w.r.t slider frame
            % Ouptut: contact   boolean value (true if there is contact)
            if abs(S_p_y) > self.slider_params.ywidth/2 %|| (S_p_x ~= -self.slider_params.xwidth/2)
                contact = false;
            else
                contact = true;
            end
            contact = true;
        end

        % Evaluate model (numerical)
        function [x_dot, mode] = eval_model(self,x,u)
            % This method can be used to evaluate the pusher-Slider non linear model x_dot = f(x,u)
            % Input: x = [x_s, y_s, theta_s (rad), r_y],
            %        u = [u_n, u_t],
            % Output: x_dot

            % State variables
            theta_s = x(3);    % rotation angle of the slider frame w.r.t. the world frame
            S_p_x = -0.034;% x(4);      % x-coordinate of the pusher position w.r.t. the slider frame S
            S_p_y = x(4);      % y-coordinate of the pusher position w.r.t. the slider frame S
            u_n = u(1);        % normal pusher velocity w.r.t. the slider frame S
            u_t = u(2);        % tangential pusher velocity w.r.t. the slider frame S

            % 2D Rotation matrix about z-axis
            R_z = helper.my_rotz(theta_s); R_z = R_z(1:2,1:2);

            % Motion Cone edges
            gamma_l = (self.slider_params.mu_sp*self.slider_params.c_ellipse^2-S_p_x*S_p_y+self.slider_params.mu_sp*S_p_x^2)/(self.slider_params.c_ellipse^2+S_p_y^2-self.slider_params.mu_sp*S_p_x*S_p_y);
            gamma_r = (-self.slider_params.mu_sp*self.slider_params.c_ellipse^2-S_p_x*S_p_y-self.slider_params.mu_sp*S_p_x^2)/(self.slider_params.c_ellipse^2+S_p_y^2+self.slider_params.mu_sp*S_p_x*S_p_y);
            v_l = [1 gamma_l]'; % left edge of the motion cone
            v_r = [1 gamma_r]'; % right edge of the motion cone


            % Evaluating u_t/u_n
            u_fract = u_t/u_n;

            % Check boundary of the motion cone
            if not(contact_surface(self,S_p_x,S_p_y))
                mode="NC";
                disp("Not in contact!");
            else
                if (u_fract <= gamma_l) && (u_fract >= gamma_r)
                    mode = "ST";
                    %                     disp("Motion Cone CHECK: Sticking Mode")
                elseif u_fract > gamma_l
                    mode = "SL";
                    %                     disp("Motion Cone CHECK: Sliding Left Mode")
                else
                    mode = "SR";
                    %                     disp("Motion Cone CHECK: Sliding Right Mode")
                end
            end


            % Model matrices
            factor_matrix = 1/(self.slider_params.c_ellipse^2+S_p_x^2+S_p_y^2);
            Q = [self.slider_params.c_ellipse^2+S_p_x^2 S_p_x*S_p_y; S_p_x*S_p_y self.slider_params.c_ellipse^2+S_p_y^2];
            %             mode = 'ST';
            switch mode
                case 'ST'
                    P = eye(2);
                    b = [-S_p_y S_p_x]';
%                                         disp('sticking mode')
                case 'SL'                   % the pusher slides on the object surface
                    P = [v_l zeros(2,1)];
                    b = [-S_p_y+gamma_l*S_p_x 0]';
%                                         disp('sliding left mode')
                case 'SR'                   % the pusher slides on the object surface
                    P = [v_r zeros(2,1)];
                    b = [-S_p_y+gamma_r*S_p_x 0]';
%                                         disp('sliding right mode')
                otherwise                   % the pusher is not in contact with the slider
                    %disp('no contact')
                    x_dot = [0 0 0 u_t]';
                    return;
            end
            c = eye(2)-factor_matrix*(Q*P+[-S_p_y; S_p_x]*b');
            F = [R_z*factor_matrix*Q*P; factor_matrix*b'; c(end,:)];
            x_dot = F*[u_n u_t]';
        end
        
        % Evaluate model with variable shape (numerical)
        function [x_dot, mode] = eval_model_variable_shape(self,x,u)
            % This method returns the symbolic expression of the nonlinear pusher-slider model x_dot = f(x,u)
            % Output: x_dot

            import casadi.*

            % named symbolic variables
            theta = x(3); % rotation angle of the slider frame w.r.t. the world frame [rad]
            s = x(4);         % curvilinear abscissa of the pusher position w.r.t. the slider frame S [m]

            u_n = u(1);     % normal pusher velocity w.r.t. slider frame S [m/s]
            u_t = u(2);     % tangential pusher velocity w.r.t. slider frame S [m/s]

            % s -> [S_p_x, S_p_y] slider frame
            if s < 0
                s = self.SP.b + s;
            end
            s = mod(s,self.SP.b);
            
            S_p = self.SP.evalSpline(self.SP.FC,s);

            % conversion [S_p_x, S_p_y] and theta to normal-tangential
            S_R_NT = full(self.SP.R_NT_fun(s));%self.SP.evalSpline(self.SP.R_NT_fun,s);
            NT_p = S_R_NT'*S_p';
            S_p_x = NT_p(1);
            S_p_y = NT_p(2);

            %             S_theta_NT = atan2(S_R_NT(2,1),S_R_NT(1,1));
            %             W_theta_NT = theta + S_theta_NT;

            % Model matrices
            sin_theta = sin(theta);
            cos_theta = cos(theta);
            c_ellipse = self.slider_params.c_ellipse;
            mu_sp = self.slider_params.mu_sp;
            factor_matrix = 1/(c_ellipse^2+S_p_x^2+S_p_y^2);

            % Motion cone edge
            gamma_l = (mu_sp*c_ellipse^2-S_p_x*S_p_y+mu_sp*S_p_x^2)/(c_ellipse^2+S_p_y^2-mu_sp*S_p_x*S_p_y);
            gamma_r = (-mu_sp*c_ellipse^2-S_p_x*S_p_y-mu_sp*S_p_x^2)/(c_ellipse^2+S_p_y^2+mu_sp*S_p_x*S_p_y);
            v_l = [1 gamma_l]';     % left edge of the motion cone
            v_r = [1 gamma_r]';     % right edge of the motion cone
            u_fract = u_t/u_n;

            % Model
            W_R_S = [cos_theta -sin_theta;sin_theta cos_theta];
            Q = [c_ellipse^2+S_p_x^2 S_p_x*S_p_y; S_p_x*S_p_y c_ellipse^2+S_p_y^2];
            % Sticking
            P_st = eye(2);
            b_st = [-S_p_y S_p_x]';
            c_st = eye(2)-factor_matrix*(Q*P_st+[-S_p_y; S_p_x]*b_st');
            F_st = [W_R_S*S_R_NT*factor_matrix*Q*P_st; factor_matrix*b_st'];
            x_dot_st = [F_st*[u_n u_t]'; 0];


            % Sliding left
            P_sl = [v_l zeros(2,1)];
            b_sl = [-S_p_y+gamma_l*S_p_x 0]';
            c_sl = eye(2)-factor_matrix*(Q*P_sl+[-S_p_y; S_p_x]*b_sl');

            % ---- s_dot evaluation
            S_p_dot_sl = S_R_NT*c_sl*[u_n u_t]';
            FC_dot = self.SP.evalSpline(self.SP.FC_dot,s);
%             s_dot_sl = (S_p_dot_sl(1)+S_p_dot_sl(2))/(FC_dot(1)+FC_dot(2));
            vp = c_sl*[u_n u_t]';
            s_dot_sl = vp-vp(1)*v_l;

            F_sl = [W_R_S*S_R_NT*factor_matrix*Q*P_sl; factor_matrix*b_sl'];
            x_dot_sl = [F_sl*[u_n u_t]'; s_dot_sl(2)];

            % Sliding right
            P_sr = [v_r zeros(2,1)];
            b_sr = [-S_p_y+gamma_r*S_p_x 0]';
            c_sr = eye(2)-factor_matrix*(Q*P_sr+[-S_p_y; S_p_x]*b_sr');

            % ---- s_dot evaluation
            S_p_dot_sr = S_R_NT*c_sr*[u_n u_t]';
%             s_dot_sr = (S_p_dot_sr(1)+S_p_dot_sr(2))/(FC_dot(1)+FC_dot(2));
%             s_dot_sr = (eye(2) - c_sr)*[u_n u_t]';
            vp = c_sr*[u_n u_t]';
            s_dot_sr = vp-vp(1)*v_r;

            F_sr = [W_R_S*S_R_NT*factor_matrix*Q*P_sr; factor_matrix*b_sr'];
            x_dot_sr = [F_sr*[u_n u_t]'; s_dot_sr(2)];

            if not(contact_surface(self,S_p_x,S_p_y))
                mode="NC";
                disp("Not in contact!");
            else
                if (u_fract <= gamma_l) && (u_fract >= gamma_r)
                    mode = "ST";
                    %                     disp("Motion Cone CHECK: Sticking Mode")
                elseif u_fract > gamma_l
                    mode = "SL";
                    %                     disp("Motion Cone CHECK: Sliding Left Mode")
                else
                    mode = "SR";
                    %                     disp("Motion Cone CHECK: Sliding Right Mode")
                end
            end

            switch mode
                case 'ST'
                    x_dot = x_dot_st;
%                     disp('sticking mode')
                case 'SL'                   % the pusher slides on the object surface
                    x_dot = x_dot_sl;
%                     disp('sliding left mode')
                case 'SR'                   % the pusher slides on the object surface
                    x_dot = x_dot_sr;
%                     disp('sliding right mode')
                otherwise                   % the pusher is not in contact with the slider
                    disp('no contact')
                    x_dot = [0 0 0 u_t]';
                    return;
            end


        end

        % Symbolic model (symbolic)
        function model = symbolic_model(self)
            % This method returns the symbolic expression of the nonlinear pusher-slider model x_dot = f(x,u)
            % Output: x_dot

            import casadi.*

            % named symbolic variables
            x = SX.sym('x');         % x-coordinate slider w.r.t. world frame [m]
            y = SX.sym('y');         % y-coordinate slider w.r.t. world frame [m]
            theta = SX.sym('theta'); % rotation angle of the slider frame w.r.t. the world frame [rad]
            S_p_x = -0.034;% SX.sym('S_p_x'); % x-coordinate of the pusher position w.r.t. the slider frame S [m]
            S_p_y = SX.sym('S_p_y'); % y-coordinate of the pusher position w.r.t. the slider frame S [m]

            u_n = SX.sym('u_n');         % normal pusher velocity w.r.t. slider frame S [m/s]
            u_t = SX.sym('u_t');         % tangential pusher velocity w.r.t. slider frame S [m/s]

            % (unnamed) symbolic variables
            %             sym_x = vertcat(x,y,theta,S_p_x,S_p_y);     % x state vector
            sym_x = vertcat(x,y,theta,S_p_y);
            sym_u = vertcat(u_n,u_t);                   % u control vector

            % Model matrices
            sin_theta = sin(theta);
            cos_theta = cos(theta);
            c_ellipse = self.slider_params.c_ellipse;
            mu_sp = self.slider_params.mu_sp;
            factor_matrix = 1/(c_ellipse^2+S_p_x^2+S_p_y^2);

            % Motion cone edge
            gamma_l = (mu_sp*c_ellipse^2-S_p_x*S_p_y+mu_sp*S_p_x^2)/(c_ellipse^2+S_p_y^2-mu_sp*S_p_x*S_p_y);
            gamma_r = (-mu_sp*c_ellipse^2-S_p_x*S_p_y-mu_sp*S_p_x^2)/(c_ellipse^2+S_p_y^2+mu_sp*S_p_x*S_p_y);
            v_l = [1 gamma_l]';     % left edge of the motion cone
            v_r = [1 gamma_r]';     % right edge of the motion cone
            u_fract = u_t/u_n;

            % Model
            R_z = [cos_theta -sin_theta;sin_theta cos_theta];
            Q = [c_ellipse^2+S_p_x^2 S_p_x*S_p_y; S_p_x*S_p_y c_ellipse^2+S_p_y^2];
            % Sticking
            P_st = eye(2);
            b_st = [-S_p_y S_p_x]';
            c_st = eye(2)-factor_matrix*(Q*P_st+[-S_p_y; S_p_x]*b_st');
            F_st = [R_z*factor_matrix*Q*P_st; factor_matrix*b_st'; c_st(end,:)];
            x_dot_st = F_st*[u_n u_t]';
            % Sliding left
            P_sl = [v_l zeros(2,1)];
            b_sl = [-S_p_y+gamma_l*S_p_x 0]';
            c_sl = eye(2)-factor_matrix*(Q*P_sl+[-S_p_y; S_p_x]*b_sl');
            F_sl = [R_z*factor_matrix*Q*P_sl; factor_matrix*b_sl'; c_sl(end,:)];
            x_dot_sl = F_sl*[u_n u_t]';
            % Sliding right
            P_sr = [v_r zeros(2,1)];
            b_sr = [-S_p_y+gamma_r*S_p_x 0]';
            c_sr = eye(2)-factor_matrix*(Q*P_sr+[-S_p_y; S_p_x]*b_sr');
            F_sr = [R_z*factor_matrix*Q*P_sr; factor_matrix*b_sr'; c_sr(end,:)];
            x_dot_sr = F_sr*[u_n u_t]';

            % f = Function('f',{sym_x,sym_u},{(u_fract>=gamma_r)*x_dot_st*(u_fract<=gamma_l) ...
            %     + (u_fract>gamma_l)*x_dot_sl...
            %     + (u_fract<gamma_r)*x_dot_sr} ...
            %     );

            % symbolic quintic function
            % %             s0 = SX.sym('s0');
            % %             sf = SX.sym('sf');
            % %             tau = SX.sym('tau');
            % %             s_tilde = Function('s_tilde',{tau},{6*tau^5-15*tau^4+10*tau^3});
            % %             s_t = Function('s_t',{s0,sf,tau},{s0 + (sf-s0)*s_tilde(tau)});

            % epsilon switching values
            % %             eps_sl = abs(0.5*gamma_l);
            % %             eps_sr = abs(0.5*gamma_r);

            % Symbolic function sticking decision
            % %             S_st = Function('S_st_fun',{sym_x,sym_u},{(u_fract>=gamma_r+eps_sr)*(u_fract<=gamma_l-eps_sl) ...
            % %                 + (u_fract<gamma_r+eps_sr)*(u_fract>gamma_r-eps_sr)*s_t(0,1,(u_fract-(gamma_r-eps_sr))/(2*eps_sr))...
            % %                 + (u_fract<gamma_l+eps_sl)*(u_fract>gamma_l-eps_sl)*s_t(1,0,(u_fract-(gamma_l-eps_sl))/(2*eps_sl))...
            % %                 } ...
            % %                 );
            % %
            % %             S_sl = Function('S_sl_fun',{sym_x,sym_u},{(u_fract<gamma_l+eps_sl)*(u_fract>gamma_l-eps_sl)*s_t(0,1,(u_fract-(gamma_l-eps_sl))/(2*eps_sl))...
            % %                 + (u_fract>gamma_l+eps_sl)} ...
            % %                 );
            % %
            % %             S_sr = Function('S_sr_fun',{sym_x,sym_u},{(u_fract<gamma_r+eps_sr)*(u_fract>gamma_r-eps_sr)*s_t(1,0,(u_fract-(gamma_r-eps_sr))/(2*eps_sr))...
            % %                 + (u_fract<gamma_r-eps_sr)} ...
            % %                 );

            % explicit dynamic function

            expr_f_expl = vertcat((u_fract>=gamma_r)*x_dot_st*(u_fract<=gamma_l) ...
                + (u_fract>gamma_l)*x_dot_sl...
                + (u_fract<gamma_r)*x_dot_sr);

            %             expr_f_expl = vertcat((S_st(sym_x,sym_u)*x_dot_st ...
            %                 + S_sl(sym_x,sym_u)*x_dot_sl...
            %                 + S_sr(sym_x,sym_u)*x_dot_sr));
            %
            %             expr_f_expl = vertcat((1*x_dot_st ...
            %                 + 0*x_dot_sl...
            %                 + 0*x_dot_sr));

            % implicit dynamic function
            %             expr_f_impl = expr_f_expl - sym_xdot;

            % Populate structure
            model.nx = self.nx;
            model.nu = self.nu;
            model.sym_x = sym_x;
            %             model.sym_xdot = sym_xdot;
            model.sym_u = sym_u;
            model.expr_f_expl = expr_f_expl;
            %             model.expr_f_impl = expr_f_impl;

            self.sym_model = model;
        end

        % Symbolic model for variable shape (symbolic)
        function model = symbolic_model_variable_shape(self)
            % This method returns the symbolic expression of the nonlinear pusher-slider model x_dot = f(x,u)
            % Output: x_dot

            import casadi.*

            % named symbolic variables
            x = SX.sym('x');         % x-coordinate slider w.r.t. world frame [m]
            y = SX.sym('y');         % y-coordinate slider w.r.t. world frame [m]
            theta = SX.sym('theta'); % rotation angle of the slider frame w.r.t. the world frame [rad]
            s = SX.sym('s');         % curvilinear abscissa of the pusher position w.r.t. the slider frame S [m]

            u_n = SX.sym('u_n');     % normal pusher velocity w.r.t. slider frame S [m/s]
            u_t = SX.sym('u_t');     % tangential pusher velocity w.r.t. slider frame S [m/s]

            % (unnamed) symbolic variables
            sym_x = vertcat(x,y,theta,s);  % x state vector
            sym_u = vertcat(u_n,u_t);      % u control vector

            % s -> [S_p_x, S_p_y] slider frame

            s = s + (s<0)*(self.SP.b);

            s = mod(s,self.SP.b);
            S_p = self.SP.FC(s);

            % conversion [S_p_x, S_p_y] and theta to normal-tangential
            S_R_NT = self.SP.R_NT_fun(s);
            NT_p = S_R_NT'*S_p';
            S_p_x = NT_p(1);
            S_p_y = NT_p(2);

            %             S_theta_NT = atan2(S_R_NT(2,1),S_R_NT(1,1));
            %             W_theta_NT = theta + S_theta_NT;

            % Model matrices
            sin_theta = sin(theta);
            cos_theta = cos(theta);
            c_ellipse = self.slider_params.c_ellipse;
            mu_sp = self.slider_params.mu_sp;
            factor_matrix = 1/(c_ellipse^2+S_p_x^2+S_p_y^2);

            % Motion cone edge
            gamma_l = (mu_sp*c_ellipse^2-S_p_x*S_p_y+mu_sp*S_p_x^2)/(c_ellipse^2+S_p_y^2-mu_sp*S_p_x*S_p_y);
            gamma_r = (-mu_sp*c_ellipse^2-S_p_x*S_p_y-mu_sp*S_p_x^2)/(c_ellipse^2+S_p_y^2+mu_sp*S_p_x*S_p_y);
            v_l = [1 gamma_l]';     % left edge of the motion cone
            v_r = [1 gamma_r]';     % right edge of the motion cone
            u_fract = u_t/u_n;

            % Model
            W_R_S = [cos_theta -sin_theta;sin_theta cos_theta];
            Q = [c_ellipse^2+S_p_x^2 S_p_x*S_p_y; S_p_x*S_p_y c_ellipse^2+S_p_y^2];
            % Sticking
            P_st = eye(2);
            b_st = [-S_p_y S_p_x]';
            F_st = [W_R_S*S_R_NT*factor_matrix*Q*P_st; factor_matrix*b_st'];
            x_dot_st = [F_st*[u_n u_t]'; 0];


            % Sliding left
            P_sl = [v_l zeros(2,1)];
            b_sl = [-S_p_y+gamma_l*S_p_x 0]';
            c_sl = eye(2)-factor_matrix*(Q*P_sl+[-S_p_y; S_p_x]*b_sl');

            % ---- s_dot evaluation
            NT_p_dot_sl = c_sl*[u_n u_t]';
            s_dot_sl = NT_p_dot_sl-NT_p_dot_sl(1)*v_l;
            
            F_sl = [W_R_S*S_R_NT*factor_matrix*Q*P_sl; factor_matrix*b_sl'];
            x_dot_sl = [F_sl*[u_n u_t]'; s_dot_sl(2)];

            % Sliding right
            P_sr = [v_r zeros(2,1)];
            b_sr = [-S_p_y+gamma_r*S_p_x 0]';
            c_sr = eye(2)-factor_matrix*(Q*P_sr+[-S_p_y; S_p_x]*b_sr');

            % ---- s_dot evaluation
            NT_p_dot_sr = c_sr*[u_n u_t]';
            s_dot_sr = NT_p_dot_sr-NT_p_dot_sr(1)*v_r;

            F_sr = [W_R_S*S_R_NT*factor_matrix*Q*P_sr; factor_matrix*b_sr'];
            x_dot_sr = [F_sr*[u_n u_t]'; s_dot_sr(2)];

            expr_f_expl = vertcat((u_fract>=gamma_r)*x_dot_st*(u_fract<=gamma_l) ...
                + (u_fract>gamma_l)*x_dot_sl...
                + (u_fract<gamma_r)*x_dot_sr);

            self.xdot_func = Function('xdot_func',{sym_x,sym_u},{expr_f_expl});

            % Populate structure
            model.nx = self.nx;
            model.nu = self.nu;
            model.sym_x = sym_x;
            %             model.sym_xdot = sym_xdot;
            model.sym_u = sym_u;
            model.expr_f_expl = expr_f_expl;
            %             model.expr_f_impl = expr_f_impl;

            self.sym_model = model;
        end
        
        % Evaluate symbolic model with numeric values
        function x_dot = evalModelVariableShape(self,x,u)
            x_dot = full(self.xdot_func(x,u));
        end
    end
end