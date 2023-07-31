clear variables
close all  
clc


%%%%%%%%%%%%%%%%%%%%%%%%%  PARAMETERS %%%%%%%%%%%%%%%%%%%%

% Spline order p
p = 3;

% Points of contact P
P = [0 0;
    0 0.5;
    0 1;
    0.5 1;
    1 1;
    1 0.5;
    1 0;
    0.5 0;
    0 0
    ];

n = length(P);

% Knots vector S = [s0 s1, ..., sm]
m = n+p+1-2*p;
a = 0; b = 4;
S_ = linspace(a,b,m);%m-2*p+1);
S = [a*ones(1,p) S_ b*ones(1,p)];
m = length(S);



%% %%%%%%%%%%%%%%%%%%%%% EVALUATION %%%%%%%%%%%%%%%%%%%%

SP = bspline_shape(S,P,p);

FC = SP.getSymbolicSpline(p);

%% PLOT
s_values = a:0.001:b;

FC_values = SP.evalSpline(FC, s_values);

plot(FC_values(:,1),FC_values(:,2),'*'), hold on, plot(P(:,1),P(:,2),'*')