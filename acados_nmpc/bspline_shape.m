classdef bspline_shape < handle
    properties
        S; % S = knot vector
        P; % P = control points
        p; % p = degree bspline
        n;
        m;
        s; % ascissa curvilinea
        C_sym; % symbolic spline
        FC; % function representing spline
        FC_dot; % function represeting dot spline
        FC_dot_dot; % function representing second derivative of the spline
        FC_angle_dot;
        n_fun; % normal versor
        t_fun; % tangential versor
        R_NT_fun; % rotation matrix for tangential and normal frame
        a;
        b;
        cj_1_vect;
        max_curvature;
    end

    methods

        function obj = bspline_shape(S,P,p)
            import casadi.*
            obj.p = p;
            obj.S = S;
            obj.P = P;
            obj.n = length(P);
            obj.m = length(S);
            obj.s = SX.sym('s');   %curvilinear abscess
            obj.FC = 0;
            obj.FC_dot = 0;
            obj.R_NT_fun = eye(2);
            obj.a = 0;
            obj.b = sum(vecnorm(diff(P)'));
        end

        function Ni_p = eval_bspline_sym(self, s, i, ord)
            % s = curvilinear abscess
            % S = knots vector
            % i = index evaluation
            % p = order bspline

            if self.S(i+ord+1)==self.S(i)
                Ni_p = 0;
                return;
            end

            if ord == 0
                Ni_p = (s<self.S(i+1))*(s>=self.S(i));
                return;
            end

            N_i_pprev = self.eval_bspline_sym(s,i,ord-1);
            N_inext_pprev = self.eval_bspline_sym(s,i+1,ord-1);

            if self.S(i+ord)==self.S(i)
                m1 = 0;
            else
                m1 = ((s-self.S(i))/(self.S(i+ord)-self.S(i)));
            end
            if self.S(i+ord+1) == self.S(i+1)
                m2 = 0;
            else
                m2 = ((self.S(i+ord+1)-s)/(self.S(i+ord+1)-self.S(i+1)));
            end

            Ni_p = m1*N_i_pprev + m2*N_inext_pprev;

        end

        function getSymbolicSpline(self,ord)
            import casadi.*
            % Get symbolic expression of spline given:
            C = 0;
            for ind = 1:self.n
                C = C + self.eval_bspline_sym(self.s,ind,ord)*self.P(ind,:);
            end
            self.C_sym = C;
            self.FC = Function('F_C',{self.s},{C});
        end

        function getSymboliSplineDot(self,ord)
            import casadi.*
            if ord < 1
                disp("ERROR: order p must be grater than 1")
            else
                self.cj_1_vect(1,:) = [0 0];
                C_dot = 0;
                for ii = 2 : self.n
                    if self.S(ii+ord)==self.S(ii)
                        cj_1 = 0;
                    else
                        cj_1 = ord*((self.P(ii,:)-self.P(ii-1,:))/(self.S(ii+ord)-self.S(ii)));
                    end
                    self.cj_1_vect(ii,:) = cj_1;
                    C_dot = C_dot + cj_1*self.eval_bspline_sym(self.s,ii,ord-1);
                end
            end

            self.FC_dot = Function('FC_dot',{self.s},{C_dot});
        end

        function getNormalTangentialVersors(self)
            import casadi.*
            t_ = self.FC_dot(self.s);
            tvers = t_/norm(t_);       %t_/norm(t_);
            nvers = -[-tvers(2) tvers(1)]; %-[-tvers(2) tvers(1)]
            R_NT = [nvers' tvers'];
            
            self.t_fun = Function('t_fun',{self.s},{tvers});
            self.n_fun = Function('n_fun',{self.s},{nvers});
            self.R_NT_fun = Function('R_NT_fun',{self.s},{R_NT});
        end

        function getSymbolicSplineDotDot(self,ord)
            import casadi.*
            if ord < 2
                disp("ERROR: order p must be grater than 2 to evaluate the second derivative")
            else
                C_dot_dot = 0;
                for ii = 3 : self.n
                    if abs(self.S(ii+ord-1)-self.S(ii))<1e-5
                        cj_2 = 0;
                    else
                        cj_2 = (ord-1)*((self.cj_1_vect(ii,:)-self.cj_1_vect(ii-1,:))/(self.S(ii+ord-1)-self.S(ii)));
                    end
%                     self.cj_1_vect(ii-1) = cj_1;
                    C_dot_dot = C_dot_dot + cj_2*self.eval_bspline_sym(self.s,ii,ord-2);
                end
            end
            self.FC_dot_dot = Function('FC_dot_dot',{self.s},{C_dot_dot});
        end
        
        function getSymbolicAngleCurvatures(self)
            import casadi.*
            t_vers_s = self.FC_dot(self.s);

            t_angle = atan2(t_vers_s(:,2),t_vers_s(:,1));

            self.FC_angle_dot = Function('t_angle_dot',{self.s},{gradient(t_angle,self.s)});
        end

        function t_angle_dot_values = getAngleCurvatures(self,s_values)
            s_values = mod(s_values, self.b);
            t_angle_dot_values = zeros(length(s_values),1);
            for j = 1:length(s_values)
                t_angle_dot_values(j) = full(self.FC_angle_dot(s_values(j)));
            end
        end

        function curvatures = getCurvatures(self,s_values)
            s_values = mod(s_values, self.b);
            curvatures = zeros(length(s_values),1);
            delta_01 = 0.011;%norm(self.P(2,:)-self.P(1,:));
            delta_0n = 0.011;%norm(self.P(end-1,:)-self.P(1,:));
            
            s1 = self.a + delta_01;
            s0 = self.a-delta_0n;

            sn = self.b + delta_01;
            sn_1 = self.b - delta_0n;

            for i = 1 : length(s_values)
                if (s_values(i) <= s1 && s_values(i) >= s0) 
                    y1 = norm(self.evalSpline(self.FC_dot_dot,s1));
                    y0 = norm(self.evalSpline(self.FC_dot_dot,s0));
                    curvatures(i) = (y1-y0)*(s_values(i)-s0)/(s1-s0) + y0;
                elseif (s_values(i) <= sn && s_values(i) >= sn_1)
                    yn = norm(self.evalSpline(self.FC_dot_dot,sn));
                    yn_1 = norm(self.evalSpline(self.FC_dot_dot,sn_1));
                    curvatures(i) = (yn-yn_1)*(s_values(i)-sn_1)/(sn-sn_1) + yn_1;
                else
                    curvatures(i) = norm(self.evalSpline(self.FC_dot_dot,s_values(i)));
                end
            end
        end
        
        function getMaxCurvature(self)
            s_values = self.a:0.001:self.b;
            curvatures = self.getCurvatures(s_values);
            self.max_curvature = max(curvatures);
        end

        function norm_curv = getNormalizedCurvature(self, s)
            curv = self.getCurvatures(s);
            norm_curv = curv/self.max_curvature;
        end

        function FC_val = evalSpline(self,F, s_values)
            s_values = mod(s_values,self.b);
            FC_val = zeros(length(s_values),2);
            % Evaluate numeric spline
            for k = 1:length(s_values)
                FC_val(k,:) = full(F(s_values(k)));
            end
        end

        function printSpline(self, a, b, step)
            s_values = a:step:b;
            FC_values = self.evalSpline(self.FC,s_values);
            t_vers_val = self.evalSpline(self.t_fun,s_values);
            n_vers_val = self.evalSpline(self.n_fun,s_values);

            quiver3(FC_values(1:1:end,1), FC_values(1:1:end,2), zeros(length(FC_values),1), t_vers_val(1:1:end,1),t_vers_val(1:1:end,2),zeros(length(FC_values),1),'r'), hold on
            quiver3(FC_values(1:1:end,1), FC_values(1:1:end,2), zeros(length(FC_values),1), n_vers_val(1:1:end,1),n_vers_val(1:1:end,2),zeros(length(FC_values),1),'b')
            plot(FC_values(:,1),FC_values(:,2),"LineWidth",2), plot(self.P(:,1),self.P(:,2),'*')
        end
    end
end