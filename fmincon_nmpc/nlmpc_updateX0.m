function options = nlmpc_updateX0(nlobj,options, x, MV)


    options.MV0 = MV;

    for i=1:nlobj.PredictionHorizon
        if i==1
            x_ = x';
        else
            x_ = options.X0(i-1,:);
        end
        options.X0(i,:) = x_ + nlobj.Ts*nlobj.Model.StateFcn(x_,options.MV0(i,:),options.Parameters{:})';
    end

end

