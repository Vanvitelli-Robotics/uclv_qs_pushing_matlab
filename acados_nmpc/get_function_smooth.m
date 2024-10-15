function get_function_smooth(p,S_p_x,S_p_y,eps_perc, u_fract_low, u_fract_upp, step_u)
    
    u_fract_vect = u_fract_low:step_u:u_fract_upp;
    
    S_st = zeros(length(u_fract_vect),1);
    S_sl = zeros(length(u_fract_vect),1);
    S_sr = zeros(length(u_fract_vect),1);
    

    for i = 1:length(u_fract_vect)

        u_fract = u_fract_vect(i);

        gamma_l = (p.slider_params.mu_sp*p.slider_params.c_ellipse^2-S_p_x*S_p_y+p.slider_params.mu_sp*S_p_x^2)/(p.slider_params.c_ellipse^2+S_p_y^2-p.slider_params.mu_sp*S_p_x*S_p_y);
        gamma_r = (-p.slider_params.mu_sp*p.slider_params.c_ellipse^2-S_p_x*S_p_y-p.slider_params.mu_sp*S_p_x^2)/(p.slider_params.c_ellipse^2+S_p_y^2+p.slider_params.mu_sp*S_p_x*S_p_y);
        
        eps_sl = abs(eps_perc * gamma_l);
        eps_sr = abs(eps_perc * gamma_r);

        % quintic function
        s_st_left = s_quintic(s_tilde(gamma_l-eps_sl,gamma_l+eps_sl,u_fract),1,0);
        s_st_right = s_quintic(s_tilde(gamma_r-eps_sr,gamma_r+eps_sr,u_fract),0,1);

        S_st(i) = (u_fract>=gamma_r+eps_sr)*(u_fract<=gamma_l-eps_sl) ...
            + (u_fract<gamma_r+eps_sr)*(u_fract>gamma_r-eps_sr)*s_st_right...
            + (u_fract<gamma_l+eps_sl)*(u_fract>gamma_l-eps_sl)*s_st_left;

        s_sl = s_quintic(s_tilde(gamma_l-eps_sl,gamma_l+eps_sl,u_fract),0,1);
        S_sl(i) = (u_fract<gamma_l+eps_sl)*(u_fract>gamma_l-eps_sl)*s_sl...
            + (u_fract>gamma_l+eps_sl);

        s_sr = s_quintic(s_tilde(gamma_r-eps_sr,gamma_r+eps_sr,u_fract),1,0);
        S_sr(i) = (u_fract<gamma_r+eps_sr)*(u_fract>gamma_r-eps_sr)*s_sr...
            + (u_fract<gamma_r-eps_sr);
    end

    figure, grid on
    plot(u_fract_vect,S_st,'*','Color','b'), title("S function"), hold on, grid on
    plot(u_fract_vect,S_sr,'*','Color','r'), 
    plot(u_fract_vect,S_sl,'*','Color','g'), 
    plot(gamma_r*ones(length(0:0.01:1)),0:0.01:1,'-.',"LineWidth",2,Color="black")
    plot(gamma_l*ones(length(0:0.01:1)),0:0.01:1,'-.',LineWidth=2, Color="black")
    legend('sticking','sliding right','sliding left')
    Sum_ = S_st+S_sl+S_sr;
    plot(u_fract_vect,Sum_,'*')
    hold off

    
    

end

function s_tilde_ = s_tilde(t0,tf,t)
    tau = (t-t0)/(tf-t0);
    s_tilde_ = 6*tau^5-15*tau^4+10*tau^3;
end

function s_quintic_ = s_quintic(s_tilde, s0, sf)
    s_quintic_ = s0 + (sf-s0)*s_tilde;
end