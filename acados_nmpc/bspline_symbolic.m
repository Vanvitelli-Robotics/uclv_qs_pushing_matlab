clear all
close all
clc

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

s_values = a:0.001:b;
C = 0;

%% Numeric evaluation


for k = 1:length(s_values)
    C = 0;
    for ind = 1:n
        C = C + eval_bspline(s_values(k),S,ind,p)*P(ind,:);
    end
    FC_values(k,:) = C; %full(F_C(s_values(k)));
end

%% Symbolic evaluation

% Script to generate spline C(s)
import casadi.*

% named symbolic variables
s = SX.sym('s');   %ascissa curvilinea

C = 0;
for ind = 1:n
    C = C + eval_bspline(s,S,ind,p)*P(ind,:);
end
F_C = Function('F_C',{s},{C});

%% PLOT

FC_values = zeros(length(s_values),2);
F_C = FC;

for k = 1:length(s_values)
    FC_values(k,:) = full(F_C(s_values(k)));
end

plot(FC_values(:,1),FC_values(:,2),'*'), hold on, plot(P(:,1),P(:,2),'*')











