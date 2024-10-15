
function [x,fval] = test_optimum_abscissa(s0_spline,c)
% options = optimoptions('fmincon','Display','off','UseParallel',false);
options = optimoptions('fmincon','Algorithm','active-set','Display','off','StepTolerance',1e-4,'OptimalityTolerance',1e-4);
[x, fval] = fmincon(@(x) optimum_abscissa(x,c),s0_spline,0,0,0,0,[],[],[],options);
end