%Rate function encoding the nonlinear equations
%of motion for the vibrating box.
%INPUTS:
%t: the current time
%V = [x;y;theta;dxdt;dydt;dthetadt]: state vector
%box_params: a struct containing the parameters that describe the system
%OUTPUTS
%dVdt = [dxdt;dydt;dthetadt;d2xdt2;d2ydt2;d2thetadt2]
function dVdt = box_rate_func(t,V,box_params)
    %#ok<INUSD>
    [ax,ay,atheta] = compute_accel(V(1),V(2),V(3),box_params);
    dVdt = [V(4:6);ax;ay;atheta];
end
