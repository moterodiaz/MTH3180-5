function results = test_box_modes()
    % Run existing tests, solver regression checks, and full modal validation.
    test_box;
    params = box_parameters();
    assert(isequal(size(params.l0_list),size(params.k_list)));
    assert(isequal(params.P_world,[-8,-2,5,9.5;-3.5,-3.5,-5,-0.5]));
    f = @(x) [x(1)^2+3*x(2);sin(x(1))-x(2)^2];
    x = [0.3;-0.2];
    [J,evals] = approximate_jacobian(f,x);
    assert(evals==5 && norm(J-[2*x(1),3;cos(x(1)),-2*x(2)],inf)<1e-8);
    [root,flag] = multi_newton_solver(@(z) [z(1)^2-2;z(2)-3], ...
        [1;1],struct('ftol',1e-12));
    assert(flag==1 && norm(root-[sqrt(2);3])<1e-10);
    results = analyze_box_modes(struct('make_video',false,'make_figures',false, ...
        'output_dir',tempname));
    assert(numel(results.modes)==3);
    for mode = 1:3
        entry = results.modes(mode);
        assert(norm(entry.initial_state(4:6))==0);
        assert(norm(entry.initial_state(1:3)-results.equilibrium(1:3)- ...
            entry.amplitude*results.mode_shapes(:,mode))<1e-14);
        assert(norm(entry.analytical(1,:)-entry.nonlinear(1,1:3))<1e-14);
        assert(entry.metrics.measured_cycles>=10);
    end
    % Smoke-test the new animator without touching the original MP4.
    entry = results.modes(results.selected_mode);
    opts = results.video_options;
    opts.output_file = fullfile(fileparts(opts.output_file),'animation_smoke_test.mp4');
    opts.figure_size = [960 540];
    opts.sample_frames = false;
    opts.fps = 10;
    info = animate_box(entry.t(1:100),entry.nonlinear(1:100,:),params,opts);
    assert(info.width==960 && info.height==540 && info.fps==10 && info.frame_count>=2);
    fprintf('All modal-analysis and animation smoke tests passed.\n');
end
