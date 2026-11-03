function simulate_box()
    %define system parameters (example values from the assignment)
    LW = 10; LH = 1; LG = 3;
    m = 1; Ic = (1/12)*(LH^2+LW^2);
    g = 1; k = 20; k_list = [.5*k,.5*k,2*k,5*k];
    l0 = 1.5*LG;
    Pbl_box = [-LW;-LH]/2;
    Pbr_box = [LW;-LH]/2;
    Ptl_box = [-LW;LH]/2;
    Ptr_box = [LW;LH]/2;
    boundary_pts = [Pbl_box,Pbr_box,Ptr_box,Ptl_box,Pbl_box];
    Pbl1_world = Pbl_box + [-LG;-LG];
    Pbl2_world = Pbl_box + [LG;-LG];
    Pbr1_world = Pbr_box + [0;-l0];
    Pbr2_world = Pbr_box + [l0;0];
    P_world = [Pbl1_world,Pbl2_world,Pbr1_world,Pbr2_world];
    P_box = [Pbl_box,Pbl_box,Pbr_box,Pbr_box];
    box_params = struct('m',m,'I',Ic,'g',g,'k_list',k_list, ...
        'l0_list',l0*ones(size(P_world,2)),'P_world',P_world, ...
        'P_box',P_box,'boundary_pts',boundary_pts);

    my_rate_func = @(t_in,V_in) box_rate_func(t_in,V_in,box_params);
    % Perturb equilibrium slightly to show coupled motion for many cycles.
    V0 = [0.6;0.8;0.25;0;0;0];
    tspan = [0,30];
    RK4.A = [0 0 0 0;1/2 0 0 0;0 1/2 0 0;0 0 1 0];
    RK4.B = [1/6;1/3;1/3;1/6];
    RK4.C = [0;1/2;1/2;1];
    [tlist,Vlist] = explicit_RK_fixed_step_integration(my_rate_func,tspan,V0,0.005,RK4);
    animate_box(tlist,Vlist,box_params);
end
