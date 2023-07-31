function Ni_p = eval_bspline(s, S, i, p)
    % s = curvilinear abscess
    % S = knots vector
    % i = index evaluation
    % p = order bspline

    if S(i+p+1)==S(i)
        Ni_p = 0;
        return;
    end
    
    if p == 0
        Ni_p = (s<S(i+1))*(s>=S(i));
        return;
    end

    N_i_pprev = eval_bspline(s,S,i,p-1);
    N_inext_pprev = eval_bspline(s,S,i+1,p-1);

    if S(i+p)==S(i)
        m1 = 0;
    else
        m1 = ((s-S(i))/(S(i+p)-S(i)));
    end
    if S(i+p+1) == S(i+1)
        m2 = 0;
    else
        m2 = ((S(i+p+1)-s)/(S(i+p+1)-S(i+1)));
    end

    Ni_p = m1*N_i_pprev + m2*N_inext_pprev;

end