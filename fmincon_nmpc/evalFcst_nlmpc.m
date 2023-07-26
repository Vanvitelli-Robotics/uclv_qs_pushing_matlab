function f = evalFcst_nlmpc(nlobj,x, MV, lastMV, varargin)

    %% process "ref" to get p-by-ny OV references (unscaled, default to 0)
if nargin >= 4
    ref = varargin{1};
else
    ref = [];
end
%% process "MVTarget" to get p-by-nmv MV references (unscaled, default to 0)
if nargin >= 5
    md = varargin{2};
else
    md = [];
end
%% process options (default to "nlmpcmoveopt")
if nargin >= 6
    options = varargin{3};
    validateattributes(options,{'nlmpcmoveopt'},{'scalar'},'nlmpcmove','"options"');
else
    options = nlmpcmoveopt;
end

%     coredata = nlmpc_getcoredata_(nlobj);

    coredata = getCoreData(nlobj);

    options = nlmpc_updateX0(nlobj,options,x, MV);

try
    [runtimedata, userdata, z0] = znlmpc_generateRuntimeData(...
        coredata, x, lastMV, ref, options.MVTarget, md, ...
        options.OutputWeights, options.MVWeights, options.MVRateWeights, options.ECRWeight, ...
        options.OutputMin, options.OutputMax, options.StateMin, options.StateMax, ...
        options.MVMin, options.MVMax, options.MVRateMin, options.MVRateMax, ...
        options.Parameters, options.X0, options.MV0, options.Slack0);
catch ME
    throw(ME)
end
%% get handles
handles = struct(...
    'hStateFcn', nlmpc__getFunctionHandle(nlobj.Model.StateFcn), ...
    'hOutputFcn', nlmpc__getFunctionHandle(nlobj.Model.OutputFcn), ...
    'hCostFcn', nlmpc__getFunctionHandle(nlobj.Optimization.CustomCostFcn), ...
    'hEqConFcn', nlmpc__getFunctionHandle(nlobj.Optimization.CustomEqConFcn), ...
    'hIneqConFcn', nlmpc__getFunctionHandle(nlobj.Optimization.CustomIneqConFcn), ...
    'hJacobianStateFcn', nlmpc__getFunctionHandle(nlobj.Jacobian.StateFcn), ...
    'hJacobianOutputFcn', nlmpc__getFunctionHandle(nlobj.Jacobian.OutputFcn), ...
    'hJacobianCostFcn', nlmpc__getFunctionHandle(nlobj.Jacobian.CustomCostFcn), ...
    'hJacobianEqConFcn', nlmpc__getFunctionHandle(nlobj.Jacobian.CustomEqConFcn), ...
    'hJacobianIneqConFcn', nlmpc__getFunctionHandle(nlobj.Jacobian.CustomIneqConFcn),...
    'hPassivityInputFcn', nlmpc__getFunctionHandle(nlobj.Passivity.InputFcn),...
    'hPassivityOutputFcn', nlmpc__getFunctionHandle(nlobj.Passivity.OutputFcn),...
    'hPassivityInputJacobianFcn', nlmpc__getFunctionHandle(nlobj.Passivity.InputJacobianFcn),...
    'hPassivityOutputJacobianFcn', nlmpc__getFunctionHandle(nlobj.Passivity.OutputJacobianFcn));

[f] = znlmpc_objfun(z0, coredata, runtimedata, userdata, handles);

end

function out = nlmpc__getFunctionHandle(in)
if isempty(in)
    out = [];
elseif isa(in,'function_handle')
    out = in;
else
    out = str2func(in);
end
end

