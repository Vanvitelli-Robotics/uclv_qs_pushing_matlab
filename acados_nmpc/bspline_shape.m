classdef bspline_shape
    properties
        S; % S = knot vector
        P; % P = control points
        p; % p = degree bspline
        n;
        m;
        s; % ascissa curvilinea
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

        function FC = getSymbolicSpline(self,ord)
            import casadi.*
            % Get symbolic expression of spline given:
            C = 0;
            for ind = 1:self.n
                C = C + self.eval_bspline_sym(self.s,ind,ord)*self.P(ind,:);
            end
            FC = Function('F_C',{self.s},{C});
        end

        function FC_dot = getSymboliSplineDot(self,ord)
            import casadi.*
            if ord < 1
                disp("ERROR: order p must be grater than 1")
            else
                C_dot = 0;
                for ii = 2 : self.n
                    if self.S(ii+ord)==self.S(ii)
                        cj_1 = 0;
                    else
                        cj_1 = ord*((self.P(ii)-self.P(ii-1))/(self.S(ii+ord)-self.S(ii)));
                    end

                    C_dot = C_dot + cj_1*self.eval_bspline_sym(self.s,ii,ord-1);
                end
            end

            FC_dot = Function('FC_dot',{self.s},{C_dot});
        end

        function FC_val = evalSpline(self,FC,s_values)
            FC_val = zeros(length(s_values),2);
            % Evaluate numeric spline
            for k = 1:length(s_values)
                FC_val(k,:) = full(FC(s_values(k)));
            end
        end
    end
end