function info = animate_box(tlist,Vlist,box_params,options)
    % Existing three-argument call retains its output name and 30 fps playback.
    % Modal options select a separate output, slow motion, and annotations.
    if nargin<4, options = struct(); end
    defaults = struct('output_file','vibrating_box.mp4','fps',30, ...
        'figure_size',[1920,1080],'playback_speed',1,'mode',[], ...
        'omega',[],'equilibrium',[],'legend_location','northwest', ...
        'visible','on','sample_frames',false,'focus_inset',false);
    fields = fieldnames(options);
    for j = 1:numel(fields)
        assert(isfield(defaults,fields{j}),'Unknown animation option: %s',fields{j});
        defaults.(fields{j}) = options.(fields{j});
    end
    options = defaults;
    validateattributes(tlist,{'numeric'},{'column','finite','increasing'});
    assert(size(Vlist,1)==numel(tlist) && size(Vlist,2)==6 && numel(tlist)>=2);
    validateattributes(options.fps,{'numeric'},{'scalar','positive','finite'});
    validateattributes(options.playback_speed,{'numeric'},{'scalar','positive','finite'});
    validateattributes(options.figure_size,{'numeric'},{'vector','numel',2,'positive','integer'});
    if ~isempty(options.mode)
        validateattributes(options.omega,{'numeric'},{'scalar','positive','finite'});
    end
    fps = options.fps;
    t_video = (tlist(1):options.playback_speed/fps:tlist(end))';
    if t_video(end) < tlist(end)-1e-10
        t_video(end+1,1) = tlist(end);
    end
    V_video = interp1(tlist,Vlist,t_video,'linear');
    fig = figure('Color','w','Units','pixels', ...
        'Position',[100 100 options.figure_size],'Visible',options.visible, ...
        'MenuBar','none','ToolBar','none','Resize','off');
    figure_cleanup = onCleanup(@() close(fig));
    ax = axes(fig); hold(ax,'on');
    if ~isempty(options.mode), ax.Position = [0.07 0.16 0.90 0.68]; end

    % Compute fixed limits once from the full simulated trajectory.
    corners = zeros(2,4,numel(t_video));
    for j = 1:numel(t_video)
        corners(:,:,j) = compute_rbt(V_video(j,1),V_video(j,2),V_video(j,3), ...
            box_params.boundary_pts(:,1:4));
    end
    all_pts = [reshape(corners,2,[]),box_params.P_world];
    span = max(all_pts,[],2)-min(all_pts,[],2);
    margin = max(span)*0.12+0.5;
    lim = [min(all_pts,[],2)-margin,max(all_pts,[],2)+margin];
    axis(ax,'equal'); xlim(ax,lim(1,:)); ylim(ax,lim(2,:));
    set(ax,'XLimMode','manual','YLimMode','manual','TickLabelInterpreter','latex', ...
        'FontSize',20,'LineWidth',1.2,'Color',[0.97 0.98 1]);
    grid(ax,'on'); box(ax,'on');
    xlabel(ax,'$x$ (-)','Interpreter','latex','FontSize',24);
    ylabel(ax,'$y$ (-)','Interpreter','latex','FontSize',24);
    title(ax,'Vibrating Box: Coupled Translation and Rotation','Interpreter','latex','FontSize',28);
    if ~isempty(options.equilibrium)
        eq = options.equilibrium;
        outline = compute_rbt(eq(1),eq(2),eq(3),box_params.boundary_pts);
        equilibrium_plot = plot(ax,outline(1,:),outline(2,:),'--', ...
            'Color',[0.56 0.60 0.65],'LineWidth',1.5,'DisplayName','Equilibrium');
    end

    spring_plots = cell(1,size(box_params.P_box,2));
    for i = 1:numel(spring_plots)
        spring_plots{i} = initialize_spring_plot(5,0.18);
    end
    box_pts = compute_rbt(V_video(1,1),V_video(1,2),V_video(1,3),box_params.boundary_pts);
    box_patch = patch(ax,box_pts(1,:),box_pts(2,:),[0.36 0.53 0.72], ...
        'EdgeColor',[0.12 0.2 0.3],'LineWidth',3,'FaceAlpha',0.85, ...
        'DisplayName','Box');
    for i = 1:numel(spring_plots)
        update_spring_plot(spring_plots{i},box_params.P_world(:,i), ...
            compute_rbt(V_video(1,1),V_video(1,2),V_video(1,3),box_params.P_box(:,i)));
    end
    handles = [box_patch,spring_plots{1}.line_plot];
    legend_names = {'Box','Springs'};
    if ~isempty(options.equilibrium)
        handles(end+1) = equilibrium_plot;
        legend_names{end+1} = 'Equilibrium outline';
    end
    orientation = 'vertical';
    if ~isempty(options.mode), orientation = 'horizontal'; end
    legend(ax,handles,legend_names, ...
        'Interpreter','latex','FontSize',18,'Location',options.legend_location, ...
        'Orientation',orientation);
    title_handle = get(ax,'Title');
    title_handle.String = sprintf('Vibrating Box: Coupled Motion, $t=%.2f$ (-)',t_video(1));
    if ~isempty(options.mode)
        title_handle.String = sprintf('Nonlinear vibrating box: mode %d, $\\omega=%.4f$ rad/(-)', ...
            options.mode,options.omega);
        timing = annotation(fig,'textbox',[0.08 0.86 0.84 0.05], ...
            'String','','Interpreter','latex','FontSize',22, ...
            'HorizontalAlignment','center','EdgeColor','none');
        note = sprintf(['Actual nonlinear trajectory; no motion amplification. ' ...
            '$f=%.4f$ cycles/(-).'],options.omega/(2*pi));
        annotation(fig,'textbox',[0.08 0.025 0.84 0.05], ...
            'String',note,'Interpreter','latex','FontSize',18, ...
            'HorizontalAlignment','center','EdgeColor','none');
    end
    if options.focus_inset
        assert(~isempty(options.equilibrium),'A fixed detail view needs an equilibrium.');
        ax.Position = [0.06 0.22 0.60 0.55];
        title_handle.FontSize = 23;
        detail = axes(fig,'Position',[0.75 0.32 0.21 0.40]);
        hold(detail,'on');
        focus = compute_rbt(eq(1),eq(2),eq(3),box_params.P_box(:,1));
        axis(detail,'equal');
        xlim(detail,focus(1)+[-0.30 0.30]);
        ylim(detail,focus(2)+[-0.30 0.30]);
        set(detail,'XLimMode','manual','YLimMode','manual', ...
            'TickLabelInterpreter','latex','FontSize',16,'Color',[0.97 0.98 1]);
        grid(detail,'on'); box(detail,'on');
        xlabel(detail,'$x$ (-)','Interpreter','latex','FontSize',20);
        ylabel(detail,'$y$ (-)','Interpreter','latex','FontSize',20);
        title(detail,{'Left attachment detail','(true coordinates; no motion gain)'}, ...
            'Interpreter','latex','FontSize',19);
        detail_box = patch(detail,box_pts(1,:),box_pts(2,:),[0.36 0.53 0.72], ...
            'EdgeColor',[0.12 0.2 0.3],'LineWidth',2,'FaceAlpha',0.60);
        plot(detail,outline(1,:),outline(2,:),'--','Color',[0.56 0.60 0.65],'LineWidth',1.5);
        detail_marker = plot(detail,focus(1),focus(2),'o','Color',[0.8 0.15 0.12], ...
            'MarkerFaceColor',[0.8 0.15 0.12],'MarkerSize',8);
        axes(ax);
    end
    vid = VideoWriter(options.output_file,'MPEG-4');
    vid.FrameRate = fps; vid.Quality = 100; open(vid);
    video_cleanup = onCleanup(@() close(vid));
    sample_indices = [];
    if options.sample_frames && ~isempty(options.mode)
        period = 2*pi/options.omega;
        sample_times = tlist(1)+[0,0.25,0.5,0.75,1]*period;
        sample_indices = unique([round((sample_times-tlist(1))*fps/options.playback_speed)+1,numel(t_video)]);
        sample_indices = sample_indices(sample_indices<=numel(t_video));
    end
    for j = 1:numel(t_video)
        pose = V_video(j,1:3);
        points_world = compute_rbt(pose(1),pose(2),pose(3),box_params.P_box);
        for i = 1:numel(spring_plots)
            update_spring_plot(spring_plots{i},box_params.P_world(:,i),points_world(:,i));
        end
        box_pts = compute_rbt(pose(1),pose(2),pose(3),box_params.boundary_pts);
        set(box_patch,'XData',box_pts(1,:),'YData',box_pts(2,:));
        if options.focus_inset
            set(detail_box,'XData',box_pts(1,:),'YData',box_pts(2,:));
            set(detail_marker,'XData',points_world(1,1),'YData',points_world(2,1));
        end
        if isempty(options.mode)
            title_handle.String = sprintf('Vibrating Box: Coupled Motion, $t=%.2f$ (-)',t_video(j));
        else
            if options.playback_speed<1
                playback = sprintf('Slow motion: $%.2f\\times$ slower',1/options.playback_speed);
            else
                playback = sprintf('Playback: $%.3f$ (-)/s',options.playback_speed);
            end
            timing.String = sprintf('$t=%.3f$ (-), cycle $%.2f$ of $%.0f$; %s', ...
                t_video(j),(t_video(j)-tlist(1))*options.omega/(2*pi), ...
                (tlist(end)-tlist(1))*options.omega/(2*pi),playback);
        end
        drawnow;
        frame = getframe(fig);
        % Retina displays may return device pixels rather than figure pixels.
        if size(frame.cdata,2)~=options.figure_size(1) || size(frame.cdata,1)~=options.figure_size(2)
            rows = round(linspace(1,size(frame.cdata,1),options.figure_size(2)));
            cols = round(linspace(1,size(frame.cdata,2),options.figure_size(1)));
            frame.cdata = frame.cdata(rows,cols,:);
        end
        writeVideo(vid,frame);
        if ismember(j,sample_indices)
            [folder,stem] = fileparts(options.output_file);
            imwrite(frame.cdata,fullfile(folder,sprintf('%s_frame_%04d.png',stem,j)));
        end
    end
    clear video_cleanup
    info = struct('file',options.output_file,'width',options.figure_size(1), ...
        'height',options.figure_size(2),'fps',fps,'frame_count',numel(t_video), ...
        'duration',numel(t_video)/fps,'simulation_duration',t_video(end)-t_video(1), ...
        'playback_speed',options.playback_speed,'axis_limits',lim, ...
        'sample_frame_indices',sample_indices);
end

%updates spring plotting object so that spring is plotted
%with ends located at points P1 and P2
function update_spring_plot(spring_plot_struct,P1,P2)
    dP = P2-P1;
    len = norm(dP);
    if len == 0, return; end
    R = [dP(1),-dP(2)/len;dP(2),dP(1)/len];
    plot_pts = R*spring_plot_struct.zig_zag;
    set(spring_plot_struct.line_plot,'xdata',plot_pts(1,:)+P1(1), ...
        'ydata',plot_pts(2,:)+P1(2));
    set(spring_plot_struct.point_plot,'xdata',[P1(1),P2(1)], ...
        'ydata',[P1(2),P2(2)]);
end

%create a struct containing plotting info for a single spring
%INPUTS:
%num_zigs: number of zig zags in spring drawing
%w: width of the spring drawing
function spring_plot_struct = initialize_spring_plot(num_zigs,w)
    zig_ending = [.25,.75,1;-1,1,0];
    zig_zag = zeros(2,3+3*num_zigs);
    zig_zag(:,1) = [-.5;0];
    zig_zag(:,end) = [num_zigs+.5;0];
    for n = 0:(num_zigs-1)
        zig_zag(:,(3+3*n):2+3*(n+1)) = zig_ending+[n,n,n;0,0,0];
    end
    zig_zag(1,:) = (zig_zag(1,:)-zig_zag(1,1))/(zig_zag(1,end)-zig_zag(1,1));
    zig_zag(2,:) = zig_zag(2,:)*w;
    spring_plot_struct.zig_zag = zig_zag;
    spring_plot_struct.line_plot = plot(0,0,'k','LineWidth',2,'DisplayName','Springs');
    spring_plot_struct.point_plot = plot(0,0,'o','Color',[0.8 0.15 0.12], ...
        'MarkerFaceColor',[0.8 0.15 0.12],'MarkerSize',7,'HandleVisibility','off');
end
