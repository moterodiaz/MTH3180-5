% Basic checks for the vibrating box functions.
F0 = compute_spring_force(10,2,[0;0],[2;0]);
assert(norm(F0)<1e-12);
Fstretch = compute_spring_force(10,2,[0;0],[3;0]);
assert(norm(Fstretch-[-10;0])<1e-12);
Fcompress = compute_spring_force(10,2,[0;0],[1;0]);
assert(norm(Fcompress-[10;0])<1e-12);
P = [1 0;-2 3];
assert(norm(compute_rbt(0,0,0,P)-P,'fro')<1e-12);
assert(norm(compute_rbt(0,0,pi/2,P)-[-P(2,:);P(1,:)],'fro')<1e-12);

% Check energy conservation for a simple vertical spring oscillator.
params = struct('m',1,'I',1,'g',0,'k_list',10,'l0_list',1, ...
    'P_world',[0;0],'P_box',[0;0]);
rate = @(t,V) box_rate_func(t,V,params);
RK4.A = [0 0 0 0;1/2 0 0 0;0 1/2 0 0;0 0 1 0];
RK4.B = [1/6;1/3;1/3;1/6]; RK4.C = [0;1/2;1/2;1];
[~,Y] = explicit_RK_fixed_step_integration(rate,[0,10],[0;1.1;0;0;0;0],0.002,RK4);
E = 0.5*(Y(:,4).^2+Y(:,5).^2) + 0.5*10*(sqrt(Y(:,1).^2+Y(:,2).^2)-1).^2;
assert((max(E)-min(E))/max(E)<1e-3);
fprintf('All vibrating box tests passed.\n');
