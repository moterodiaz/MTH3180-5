function [t_list,X_list,h_avg,num_evals] = explicit_RK_fixed_step_integration(rate_func_in,tspan,X0,h_ref,BT_struct)
    n_steps = max(1,ceil((tspan(2)-tspan(1))/h_ref));
    t_list = linspace(tspan(1),tspan(2),n_steps+1)';
    h_avg = (tspan(2)-tspan(1))/n_steps;
    X_list = zeros(n_steps+1,numel(X0));
    X_list(1,:) = X0(:)';
    num_evals = 0;
    for n = 1:n_steps
        [XB,evals] = explicit_RK_step(rate_func_in,t_list(n),X_list(n,:)',h_avg,BT_struct);
        X_list(n+1,:) = XB';
        num_evals = num_evals + evals;
    end
end

function [XB,num_evals] = explicit_RK_step(rate_func_in,t,XA,h,BT_struct)
    XA = XA(:);
    s = length(BT_struct.C);
    K = zeros(numel(XA),s);
    for i = 1:s
        X_stage = XA + h*K*BT_struct.A(i,:)';
        K(:,i) = rate_func_in(t+BT_struct.C(i)*h,X_stage);
    end
    XB = XA + h*K*BT_struct.B(:);
    num_evals = s;
end
