clear variables
close all
clc

% NMPC Parameters
nx = 5; % number of state variables
ny = 5; % number of output variables
nu = 2; % number of control variables
Ts = 0.05; % sample time
Hp = 10; % prediction horizon   
Hu = Hp; % control horizon

% Create nlmpc object
nlobj = nlmpc(nx,ny,nu);
nlobj.Ts = Ts;
nlobj.PredictionHorizon = Hp;
nlobj.ControlHorizon = Hu;
nlobj.Model.StateFcn = "";
nlobj.Model.IsContinousTime = true;
nlobj.Model.OutputFcn = "";
% nlobj.Jacobian.OutputFcn = ;

% Define Cost and Constraints 
nlobj.Optimization.ReplaceStandardCost = false;
