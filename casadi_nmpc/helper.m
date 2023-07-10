classdef helper
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
   end

end