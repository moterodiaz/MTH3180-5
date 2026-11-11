function checks = verify_mode_artifacts(output_dir)
    % Check saved figures and decode the completed MP4 using MATLAB.
    if nargin<1
        output_dir = fileparts(mfilename('fullpath'));
    end
    data = load(fullfile(output_dir,'modal_results.mat'),'results');
    results = data.results;
    checks.png_sizes = zeros(3,2);
    for mode = 1:3
        stem = fullfile(output_dir,sprintf('mode_%d_comparison',mode));
        info = imfinfo([stem,'.png']);
        assert(info.Width>=3000 && info.Height>=2000,'Comparison PNG resolution is too low.');
        checks.png_sizes(mode,:) = [info.Width,info.Height];
        fig = openfig([stem,'.fig'],'invisible');
        cleanup = onCleanup(@() close(fig));
        axes_list = findall(fig,'Type','axes');
        assert(numel(axes_list)==3);
        for ax = axes_list'
            assert(strcmp(ax.TickLabelInterpreter,'latex'));
            assert(strcmp(ax.YLabel.Interpreter,'latex'));
            assert(numel(findall(ax,'Type','line'))==2);
            assert(strcmp(ax.XLimMode,'manual') && strcmp(ax.YLimMode,'manual'));
        end
        lgd = findall(fig,'Type','legend');
        assert(isscalar(lgd) && strcmp(lgd.Interpreter,'latex'));
        assert(strcmp(lgd.Layout.Tile,'south'),'Shared legend is not outside the plots.');
        clear cleanup
    end
    reader = VideoReader(results.video_file);
    checks.width = reader.Width;
    checks.height = reader.Height;
    checks.fps = reader.FrameRate;
    checks.duration = reader.Duration;
    checks.frame_count = reader.NumFrames;
    assert(checks.width==1920 && checks.height==1080 && abs(checks.fps-60)<1e-9);
    assert(checks.frame_count==results.video.frame_count);
    assert(abs(checks.duration-checks.frame_count/checks.fps)<1/checks.fps);
    checks.measured_cycles = results.modes(results.selected_mode).metrics.measured_cycles;
    assert(checks.measured_cycles>=10);
    sample_indices = results.video.sample_frame_indices;
    for index = sample_indices
        frame = read(reader,index);
        assert(isequal(size(frame),[1080,1920,3]));
        imwrite(frame,fullfile(output_dir,sprintf('decoded_frame_%04d.png',index)));
    end
    checks.decoded_frame_indices = sample_indices;
    fid = fopen(fullfile(output_dir,'artifact_validation.json'),'w');
    assert(fid~=-1);
    cleanup = onCleanup(@() fclose(fid));
    fprintf(fid,'%s\n',jsonencode(checks,'PrettyPrint',true));
    disp(checks);
    fprintf('All saved figure and video checks passed.\n');
end
