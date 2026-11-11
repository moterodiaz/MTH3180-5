function simulate_box()
    box_params = box_parameters();

    my_rate_func = @(t_in,V_in) box_rate_func(t_in,V_in,box_params);
    % Perturb equilibrium slightly to show coupled motion for many cycles.
    V0 = [0.6;0.8;0.25;0;0;0];
    tspan = [0,30];
    RK4 = rk4_tableau();
    [tlist,Vlist] = explicit_RK_fixed_step_integration(my_rate_func,tspan,V0,0.005,RK4);
    animate_box(tlist,Vlist,box_params);
end
