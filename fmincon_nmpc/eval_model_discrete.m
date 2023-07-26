function xk = eval_model_discrete(x_dot, x0, Ts)
    xk = x0 + Ts*x_dot;
end

