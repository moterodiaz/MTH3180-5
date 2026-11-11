function RK4 = rk4_tableau()
    % Classical fourth-order explicit Runge-Kutta method.
    RK4.A = [0 0 0 0;1/2 0 0 0;0 1/2 0 0;0 0 1 0];
    RK4.B = [1/6;1/3;1/3;1/6];
    RK4.C = [0;1/2;1/2;1];
end
