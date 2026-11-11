function params = box_parameters()
    % Assignment 5 example geometry, shared by both simulation drivers.
    LW = 10; LH = 1; LG = 3;
    m = 1; Ic = (1/12)*(LH^2+LW^2);
    g = 1; k = 20; k_list = [.5*k,.5*k,2*k,5*k];
    l0 = 1.5*LG;
    Pbl_box = [-LW;-LH]/2;
    Pbr_box = [LW;-LH]/2;
    Ptl_box = [-LW;LH]/2;
    Ptr_box = [LW;LH]/2;
    boundary_pts = [Pbl_box,Pbr_box,Ptr_box,Ptl_box,Pbl_box];
    P_world = [Pbl_box+[-LG;-LG],Pbl_box+[LG;-LG], ...
        Pbr_box+[0;-l0],Pbr_box+[l0;0]];
    P_box = [Pbl_box,Pbl_box,Pbr_box,Pbr_box];
    params = struct('m',m,'I',Ic,'g',g,'k_list',k_list, ...
        'l0_list',l0*ones(1,size(P_world,2)),'P_world',P_world, ...
        'P_box',P_box,'boundary_pts',boundary_pts);
end
