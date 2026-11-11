
clear; clc; close all;

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

%define system parameters
box_params = struct();
box_params.m = m;
box_params.I = Ic;
box_params.g = g;
box_params.k_list = k_list;
box_params.l0_list = l0*ones(size(P_world,2));
box_params.P_world = P_world;
box_params.P_box = P_box;
box_params.boundary_pts = boundary_pts;

my_rate_func = @(t_in,V_in) box_rate_func(t_in,V_in,box_params);

Veq = zeros(6,1);
for n = 1:20
    F = my_rate_func(0,Veq);
    J = zeros(3,3);
    for j = 1:3
        h = zeros(6,1);
        h(j) = 1e-6;
        Fp = my_rate_func(0,Veq+h);
        Fm = my_rate_func(0,Veq-h);
        J(:,j) = (Fp(4:6)-Fm(4:6))/(2e-6);
    end
    if norm(F(4:6)) < 1e-10
        break
    end
    Veq(1:3) = Veq(1:3)-J\F(4:6);
end

if norm(my_rate_func(0,Veq)) > 1e-7
    error('Equilibrium calculation did not converge')
end

J_approx = zeros(6,6);
for j = 1:6
    h = zeros(6,1);
    h(j) = 1e-6;
    J_approx(:,j) = (my_rate_func(0,Veq+h)-my_rate_func(0,Veq-h))/(2e-6);
end

V_eq = Veq;
my_linear_rate = @(t_in,V_in) J_approx*(V_in-V_eq);

dx0 = 0.5;
dy0 = 0.3;
dtheta0 = 0.07;
vx0 = 0;
vy0 = 0;
vtheta0 = 0;

epsilon_list = [0.05,1.5];
case_names = {'Near equilibrium','Far from equilibrium'};
file_names = {'box_near_equilibrium.png','box_far_equilibrium.png'};
tspan = 0:0.002:4;
options = odeset('RelTol',1e-9,'AbsTol',1e-11);

for j = 1:length(epsilon_list)
    epsilon = epsilon_list(j);
    V0 = Veq + epsilon*[dx0;dy0;dtheta0;vx0;vy0;vtheta0];

    [tlist_nonlinear,Vlist_nonlinear] = ode45(my_rate_func,tspan,V0,options);
    [tlist_linear,Vlist_linear] = ode45(my_linear_rate,tspan,V0,options);

    fig = figure('Color','w','Position',[100 100 1100 850]);
    tl = tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
    ylabels = {'$x$ (m)','$y$ (m)','$\theta$ (rad)'};
    plot_titles = {'Horizontal position','Vertical position','Orientation'};
    ax = gobjects(3,1);

    for i = 1:3
        ax(i) = nexttile;
        p1 = plot(tlist_nonlinear,Vlist_nonlinear(:,i),'b-','LineWidth',1.8);
        hold on
        p2 = plot(tlist_linear,Vlist_linear(:,i),'r--','LineWidth',1.8);
        if i == 1
            legend_lines = [p1 p2];
        end
        title(plot_titles{i},'Interpreter','latex','FontSize',16);
        xlabel('$t$ (s)','Interpreter','latex','FontSize',15);
        ylabel(ylabels{i},'Interpreter','latex','FontSize',15);
        set(gca,'TickLabelInterpreter','latex','FontSize',13);
        xlim([0 4]);
        grid on
    end

    lgd = legend(ax(1),legend_lines,{'Nonlinear','Linear'}, ...
        'Interpreter','latex','FontSize',14,'Orientation','horizontal');
    lgd.Layout.Tile = 'south';
    sgtitle(tl,sprintf('%s ($\\epsilon = %.2f$)',case_names{j},epsilon), ...
        'Interpreter','latex','FontSize',18);
    exportgraphics(fig,file_names{j},'Resolution',300);
end
