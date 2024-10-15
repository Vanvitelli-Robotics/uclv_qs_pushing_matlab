function slider = object_selection(obj)
    switch obj
        case 'santal'
            slider.mu_sg = 0.32;                                  % friction coefficient between slider and ground
            slider.mu_sp = 0.19;                                 % friction coefficient between slider and pusher
            slider.xwidth = 0.068;                               % width of the slider along x-direction [m]
            slider.ywidth = 0.082;                                % width of the slider along y-direction [m]
            slider.area = slider.xwidth * slider.ywidth;          % slider area [m^2]
            slider.m = 0.2875;                                    % slider mass [kg]
            slider.tau_max = 0.0251;
            slider.cad_model_path = "cad_santal_centered_scaled_rotated_reduced.stl";
            slider.pcl_path = "planar_surface_santal_36_uniformed.ply";%
        case 'balea'
            slider.mu_sg = 0.35;
            slider.mu_sp = 0.20;
            slider.xwidth = 0.071;                               % width of the slider along x-direction [m]
            slider.ywidth = 0.071;                                % width of the slider along y-direction [m]
            slider.area = slider.xwidth * slider.ywidth;          % slider area [m^2]
            slider.m = 0.1713;
            slider.tau_max = 0.0042;
            slider.cad_model_path = "Balea_cad_model v1.stl";
            slider.pcl_path = 'Balea_cad_model_planar_surface_36.ply';
        case 'montana'
            slider.mu_sg = 0.20;
            slider.mu_sp = 0.10;
            slider.xwidth = 0.057;                               % width of the slider along x-direction [m]
            slider.ywidth = 0.101;                                % width of the slider along y-direction [m]
            slider.area = slider.xwidth * slider.ywidth;          % slider area [m^2]
            slider.m = 0.2467;
            slider.tau_max = 0.0101;
            slider.cad_model_path = "Montana_cad_model.stl";
            slider.pcl_path = 'Montana_cad_model_planar_section_34.ply';
        case 'pulirapid'
            slider.mu_sg = 0.22;                                  % friction coefficient between slider and ground
            slider.mu_sp = 0.1;                                 % friction coefficient between slider and pusher
            slider.xwidth = 0.13;                               % width of the slider along x-direction [m]
            slider.ywidth = 0.23;                                % width of the slider along y-direction [m]
            slider.area = slider.xwidth * slider.ywidth;          % slider area [m^2]
            slider.m = 0.500;                                    % slider mass [kg]
            slider.tau_max = 0.0251;
            slider.cad_model_path = "pulirapid_ricarica_simplified.stl";
            slider.pcl_path = "pulirapid_ricarica_test_curvatura2_ply.ply";%
        otherwise
            disp("Invalid object! Please, chose between: santal, balea, montana")
            return;
    end
end