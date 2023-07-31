classdef bspline_shape
    properties
        S; % S = knot vector
        P; % P = control points
        p; % p = degree bspline
        n;
        m;
    end

    methods

        function obj = bspline_shape(S,P,p)
            obj.p = p;
            obj.S = S;
            obj.P = P;
            obj.n = length(P);
            obj.m = length(S);
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
            % Get symbolic expression of spline given:

            import casadi.*
            s = SX.sym('s');   %curvilinear abscess
            C = 0;
            for ind = 1:self.n
                C = C + self.eval_bspline_sym(s,ind,ord)*self.P(ind,:);
            end
            FC = Function('F_C',{s},{C});
        end

        %         function FC_dot = getSymboliSplineDot(P,S,p)
        %             cj = p*(()/());
        %         end

        function FC_val = evalSpline(self,FC,s_values)
            FC_val = zeros(length(s_values),2);
            % Evaluate numeric spline
            for k = 1:length(s_values)
                FC_val(k,:) = full(FC(s_values(k)));
            end
        end
    end
end