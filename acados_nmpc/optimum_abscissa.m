function f = optimum_abscissa(x,c)
    f = norm([c.S_p_x c.S_p_y] - full(c.fc_spline(x)))^2;
end