function energy = box_energy(V,params)
    % Total mechanical energy for rows [x y theta vx vy vtheta].
    energy = 0.5*params.m*sum(V(:,4:5).^2,2) + ...
        0.5*params.I*V(:,6).^2 + params.m*params.g*V(:,2);
    for j = 1:size(V,1)
        points = compute_rbt(V(j,1),V(j,2),V(j,3),params.P_box);
        lengths = vecnorm(points-params.P_world);
        energy(j) = energy(j)+0.5*sum(params.k_list.*(lengths-params.l0_list).^2);
    end
end
