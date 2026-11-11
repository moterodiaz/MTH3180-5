function results = analyze_box_modes(options)
    % Reproduce modal analysis, validation, report figures, and optional video.
    % analyze_box_modes(struct('make_video',false)) skips only video rendering.
    if nargin < 1, options = struct(); end
    defaults = struct('output_dir',fileparts(mfilename('fullpath')), ...
        'make_figures',true,'make_video',true);
    names = fieldnames(options);
    for j = 1:numel(names)
        assert(isfield(defaults,names{j}),'Unknown analysis option: %s',names{j});
        defaults.(names{j}) = options.(names{j});
    end
    options = defaults;
    if ~isfolder(options.output_dir), mkdir(options.output_dir); end
    params = box_parameters();
    rate = @(t,V) box_rate_func(t,V,params);
    root_rate = @(V) rate(0,V);
    [Veq,flag,root_evals] = multi_newton_solver(root_rate,zeros(6,1), ...
        struct('ftol',1e-12,'dxmin',1e-14,'max_iter',100));
    assert(flag==1 && norm(root_rate(Veq),inf)<1e-10,'Equilibrium solve failed.');
    [A,jac_evals] = approximate_jacobian(root_rate,Veq);
    Q = -A(4:6,1:3);
    M = diag([params.m,params.m,params.I]);
    [Phi,D] = eig(Q);
    assert(max(abs(imag(diag(D))))<1e-9 && min(real(diag(D)))>0, ...
        'Equilibrium is not a stable oscillatory equilibrium.');
    [lambda,order] = sort(real(diag(D)));
    Phi = real(Phi(:,order));
    for j = 1:3
        Phi(:,j) = Phi(:,j)/sqrt(Phi(:,j)'*M*Phi(:,j));
        [~,pivot] = max(abs(Phi(:,j)));
        if Phi(pivot,j)<0, Phi(:,j) = -Phi(:,j); end
    end
    omega = sqrt(lambda);
    periods = 2*pi./omega;
    jac_reference = five_point_jacobian(root_rate,Veq,1e-4);
    validation = struct('equilibrium_residual',norm(root_rate(Veq),inf), ...
        'jacobian_structure_error',norm(A-[zeros(3),eye(3);-Q,zeros(3)],inf), ...
        'jacobian_reference_relative_error',norm(A-jac_reference,'fro')/norm(A,'fro'), ...
        'mass_stiffness_symmetry_error',norm(M*Q-Q'*M,'fro')/norm(M*Q,'fro'), ...
        'mass_orthogonality_error',norm(Phi'*M*Phi-eye(3),'fro'), ...
        'eigenpair_relative_residual',norm(Q*Phi-Phi*diag(lambda),'fro')/norm(Q*Phi,'fro'));
    assert(validation.jacobian_structure_error<1e-10);
    assert(validation.jacobian_reference_relative_error<1e-7);
    assert(validation.mass_stiffness_symmetry_error<1e-7);
    assert(validation.mass_orthogonality_error<1e-7);
    assert(validation.eigenpair_relative_residual<1e-12);
    RK4 = rk4_tableau();
    num_cycles = 12;
    [~,rest] = explicit_RK_fixed_step_integration(rate,[0,num_cycles*periods(1)], ...
        Veq,periods(3)/160,RK4);
    validation.equilibrium_simulation_max_drift = max(abs(rest-Veq'),[],'all');
    assert(validation.equilibrium_simulation_max_drift<1e-9);
    results = struct('parameters',params,'equilibrium',Veq,'A',A,'Q',Q, ...
        'mass_matrix',M,'mode_shapes',Phi,'omega',omega, ...
        'frequency',omega/(2*pi),'periods',periods,'num_cycles',num_cycles, ...
        'newton_function_evaluations',root_evals,'jacobian_function_evaluations',jac_evals, ...
        'validation',validation,'matlab_version',version);
    % Search largest-to-smallest actual corner displacement, not display gain.
    targets = [0.30,0.20,0.14,0.10,0.07,0.05,0.045,0.04,0.035,0.025,0.018,0.012,0.008,0.004];
    for mode = 1:3
        duration = num_cycles*periods(mode);
        coarse_steps = ceil(duration/(periods(3)/80));
        h = duration/(2*coarse_steps);
        trial_history = [];
        found = false;
        for target = targets
            amplitude = target/corner_mode_scale(Phi(:,mode),Veq,params);
            V0 = Veq + [amplitude*Phi(:,mode);zeros(3,1)];
            [t,V] = explicit_RK_fixed_step_integration(rate,[0,duration],V0,h,RK4);
            metrics = modal_metrics(t,V,Veq,Phi,M,omega,mode,amplitude,params);
            trial_history(end+1,:) = [target,amplitude,metrics.rms_error, ...
                metrics.peak_error,metrics.leakage]; %#ok<AGROW>
            if metrics.rms_error<0.01 && metrics.peak_error<0.01 && metrics.leakage<0.01
                found = true;
                break
            end
        end
        assert(found,'No amplitude met the agreement and leakage thresholds.');
        [~,Vc] = explicit_RK_fixed_step_integration(rate,[0,duration],V0,2*h,RK4);
        [tf,Vf] = explicit_RK_fixed_step_integration(rate,[0,duration],V0,h/2,RK4);
        scales = [sqrt(diag(M));sqrt(diag(M))/omega(mode)]';
        coarse_difference = norm((Vc-V(1:2:end,:)).*scales,'fro')/sqrt(size(Vc,1));
        fine_difference = norm((V-Vf(1:2:end,:)).*scales,'fro')/sqrt(size(V,1));
        convergence_order = log2(coarse_difference/fine_difference);
        fprintf('Convergence mode %d: coarse %.9g, fine %.9g, order %.6f\n', ...
            mode,coarse_difference,fine_difference,convergence_order);
        convergence_relative_error = fine_difference/ ...
            (norm((Vf(1:2:end,:)-Veq').*scales,'fro')/sqrt(size(V,1)));
        assert(convergence_order>3.7 && convergence_order<4.3, ...
            'RK4 step-halving did not show fourth-order convergence.');
        assert(convergence_relative_error<1e-6,'RK4 error is too large.');
        metrics = modal_metrics(tf,Vf,Veq,Phi,M,omega,mode,amplitude,params);
        assert(metrics.rms_error<0.01 && metrics.peak_error<0.01 && metrics.leakage<0.01);
        assert(metrics.relative_energy_drift<1e-6);
        assert(metrics.measured_cycles>=10,'Fewer than ten measured complete cycles.');
        analytical = Veq(1:3)' + amplitude*cos(omega(mode)*tf)*Phi(:,mode)';
        entry = struct('mode',mode,'t',tf,'nonlinear',Vf,'analytical',analytical, ...
            'amplitude',amplitude,'initial_state',V0,'metrics',metrics, ...
            'trial_history',trial_history,'step_sizes',[2*h,h,h/2], ...
            'convergence_differences',[coarse_difference,fine_difference], ...
            'convergence_order',convergence_order, ...
            'convergence_relative_error',convergence_relative_error);
        results.modes(mode) = entry;
        fprintf('Mode %d: omega=%.9f, f=%.9f, amplitude=%.6g, RMS=%.5f%%, peak=%.5f%%, leakage=%.5f%%, RK order=%.4f\n', ...
            mode,omega(mode),omega(mode)/(2*pi),amplitude,100*metrics.rms_error, ...
            100*metrics.peak_error,100*metrics.leakage,convergence_order);
        if options.make_figures
            export_mode_comparison(entry,results,options.output_dir);
        end
    end
    % Largest validated true corner travel is clearest at the fixed scene scale.
    [~,selected] = max(arrayfun(@(m) m.metrics.corner_peak_to_peak,results.modes));
    results.selected_mode = selected;
    results.video_file = fullfile(options.output_dir,sprintf('mode_%d_nonlinear.mp4',selected));
    chosen = results.modes(selected);
    % A period lasts two playback seconds. The trajectory itself is unchanged.
    video_options = struct('output_file',results.video_file,'fps',60, ...
        'figure_size',[1920,1080],'playback_speed',periods(selected)/2, ...
        'mode',selected,'omega',omega(selected),'equilibrium',Veq, ...
        'legend_location','southoutside','visible','off','sample_frames',true,'focus_inset',true);
    results.video_options = video_options;
    results.video_rendered = options.make_video;
    save(fullfile(options.output_dir,'modal_results.mat'),'results');
    write_modal_validation(results,options.output_dir);
    if options.make_video
        results.video = animate_box(chosen.t,chosen.nonlinear,params,video_options);
        save(fullfile(options.output_dir,'modal_results.mat'),'results');
        write_modal_validation(results,options.output_dir);
    end
end

function scale = corner_mode_scale(phi,Veq,params)
    theta = Veq(3);
    dR = [-sin(theta),-cos(theta);cos(theta),-sin(theta)];
    velocities = phi(1:2)+phi(3)*dR*params.boundary_pts(:,1:4);
    scale = max(vecnorm(velocities));
end

function J = five_point_jacobian(fun,x,h)
    % Independent higher-order reference, not the reused Newton differentiator.
    J = zeros(numel(fun(x)),numel(x));
    for j = 1:numel(x)
        delta = zeros(size(x)); delta(j) = h;
        J(:,j) = (-fun(x+2*delta)+8*fun(x+delta)-8*fun(x-delta)+fun(x-2*delta))/(12*h);
    end
end

function metrics = modal_metrics(t,V,Veq,Phi,M,omega,mode,amplitude,params)
    displacement = V(:,1:3)-Veq(1:3)';
    prediction = amplitude*cos(omega(mode)*t)*Phi(:,mode)';
    weighted_error = (displacement-prediction).*sqrt(diag(M))';
    weighted_prediction = prediction.*sqrt(diag(M))';
    modal_position = (Phi\displacement')';
    modal_velocity = (Phi\V(:,4:6)')';
    phase_amplitudes = sqrt(modal_position.^2+(modal_velocity./omega').^2);
    other = setdiff(1:3,mode);
    energy = box_energy(V,params);
    crossings = find(modal_position(1:end-1,mode)<0 & modal_position(2:end,mode)>=0);
    crossing_times = t(crossings)-modal_position(crossings,mode).* ...
        (t(crossings+1)-t(crossings))./ ...
        (modal_position(crossings+1,mode)-modal_position(crossings,mode));
    assert(numel(crossing_times)>=2,'Insufficient cycles for frequency measurement.');
    measured_omega = 2*pi/mean(diff(crossing_times));
    base = compute_rbt(Veq(1),Veq(2),Veq(3),params.boundary_pts(:,1:4));
    corner_displacement = zeros(size(V,1),8);
    for j = 1:size(V,1)
        corners = compute_rbt(V(j,1),V(j,2),V(j,3),params.boundary_pts(:,1:4));
        corner_displacement(j,:) = reshape(corners-base,1,[]);
    end
    corner_ranges = reshape(max(corner_displacement)-min(corner_displacement),2,[]);
    metrics = struct('rms_error',norm(weighted_error,'fro')/norm(weighted_prediction,'fro'), ...
        'peak_error',max(vecnorm(weighted_error,2,2))/abs(amplitude), ...
        'coordinate_rms_error',sqrt(mean((displacement-prediction).^2))./ ...
            max(abs(amplitude*Phi(:,mode)'),1e-12), ...
        'leakage',max(vecnorm(phase_amplitudes(:,other),2,2))/abs(amplitude), ...
        'relative_energy_drift',(max(energy)-min(energy))/ ...
            (0.5*omega(mode)^2*amplitude^2), ...
        'measured_omega',measured_omega, ...
        'relative_frequency_error',abs(measured_omega/omega(mode)-1), ...
        'measured_cycles',floor((t(end)-t(1))*measured_omega/(2*pi)), ...
        'crossing_times',crossing_times,'corner_peak_to_peak',max(vecnorm(corner_ranges)));
end
