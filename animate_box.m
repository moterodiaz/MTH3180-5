function animate_box(tlist,Vlist,box_params)
    fps = 30;
    t_video = (tlist(1):1/fps:tlist(end))';
    if t_video(end) < tlist(end)-1e-10
        t_video(end+1,1) = tlist(end);
    end
    V_video = interp1(tlist,Vlist,t_video,'linear');
    fig = figure('Color','w','Units','pixels','Position',[100 100 1920 1080]);
    ax = axes(fig); hold(ax,'on');

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
    legend(ax,[box_patch,spring_plots{1}.line_plot],{'Box','Springs'}, ...
        'Interpreter','latex','FontSize',18,'Location','northwest');
    title_handle = get(ax,'Title');
    title_handle.String = sprintf('Vibrating Box: Coupled Motion, $t=%.2f$ (-)',t_video(1));
    vid = VideoWriter('vibrating_box.mp4','MPEG-4');
    vid.FrameRate = fps; vid.Quality = 100; open(vid);
    for j = 1:numel(t_video)
        pose = V_video(j,1:3);
        points_world = compute_rbt(pose(1),pose(2),pose(3),box_params.P_box);
        for i = 1:numel(spring_plots)
            update_spring_plot(spring_plots{i},box_params.P_world(:,i),points_world(:,i));
        end
        box_pts = compute_rbt(pose(1),pose(2),pose(3),box_params.boundary_pts);
        set(box_patch,'XData',box_pts(1,:),'YData',box_pts(2,:));
        title_handle.String = sprintf('Vibrating Box: Coupled Motion, $t=%.2f$ (-)',t_video(j));
        drawnow;
        writeVideo(vid,getframe(fig));
    end
    close(vid);
    close(fig);
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
