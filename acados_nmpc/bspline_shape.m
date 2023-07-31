classdef bspline_shape
    properties
    end

    methods(Static)
        function FC = getSymbolicSpline(P,S,p)
            % P = control points
            % S = knot vector
            % p = degree bspline
            import casadi.*
            s = SX.sym('s');   %ascissa curvilinea
            C = 0;
            for ind = 1:n
                C = C + eval_bspline(s,S,ind,p)*P(ind,:);
            end
            FC = Function('F_C',{s},{C});
        end

        function FC_dot = getSymboliSplineDot()
        end

        function rp = evalSpline(FC,k)
            rp = full(FC(s_values(k)));
        end
    end
end