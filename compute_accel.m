%Computes the linear and angular acceleration of the box
%given its current position and orientation
%INPUTS:
%x: current x position of the box
%y: current y position of the box
%theta: current orientation of the box
%box_params: a struct containing the parameters that describe the system
%OUTPUTS
%ax: x acceleration of the box
%ay: y acceleration of the box
%atheta: angular acceleration of the box
function [ax,ay,atheta] = compute_accel(x,y,theta,box_params)
    P_box_world = compute_rbt(x,y,theta,box_params.P_box);
    F_total = [0;-box_params.m*box_params.g];
    torque_total = 0;
    for i = 1:numel(box_params.k_list)
        F = compute_spring_force(box_params.k_list(i),box_params.l0_list(i), ...
            box_params.P_world(:,i),P_box_world(:,i));
        F_total = F_total + F;
        r = P_box_world(:,i)-[x;y];
        torque_total = torque_total + r(1)*F(2)-r(2)*F(1);
    end
    ax = F_total(1)/box_params.m;
    ay = F_total(2)/box_params.m;
    atheta = torque_total/box_params.I;
end
