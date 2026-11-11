function export_mode_comparison(entry,results,output_dir)
    % Plot absolute generalized coordinates; use narrow limits about equilibrium.
    mode = entry.mode;
    fig = figure('Color','w','Units','pixels','Position',[50 50 1600 1200], ...
        'Visible','off');
    cleanup = onCleanup(@() close(fig));
    layout = tiledlayout(fig,3,1,'TileSpacing','compact','Padding','loose');
    title(layout,sprintf('Mode %d: nonlinear simulation and analytical prediction',mode), ...
        'Interpreter','latex','FontSize',25);
    subtitle(layout,sprintf(['$\\omega_%d=%.5f$ rad/(-), $f_%d=%.5f$ cycles/(-), ' ...
        '$\\epsilon=%.5f$; 12 predicted cycles\n' ...
        'Mass-weighted RMS error: %.3f\\%%; peak unwanted-mode leakage: %.3f\\%%'], ...
        mode,results.omega(mode),mode,results.frequency(mode),entry.amplitude, ...
        100*entry.metrics.rms_error,100*entry.metrics.leakage), ...
        'Interpreter','latex','FontSize',17);
    labels = {'$x$ (-)','$y$ (-)','$\theta$ (rad)'};
    axes_list = gobjects(3,1);
    for coordinate = 1:3
        ax = nexttile(layout); axes_list(coordinate) = ax; hold(ax,'on');
        h1 = plot(ax,entry.t,entry.analytical(:,coordinate),'-', ...
            'Color',[0.12 0.32 0.58],'LineWidth',2.8);
        h2 = plot(ax,entry.t,entry.nonlinear(:,coordinate),'--', ...
            'Color',[0.85 0.30 0.09],'LineWidth',1.8);
        set(ax,'TickLabelInterpreter','latex','FontSize',18,'LineWidth',1, ...
            'XLim',[entry.t(1),entry.t(end)],'XTick',linspace(0,entry.t(end),7));
        ax.YAxis.Exponent = 0;
        xtickformat(ax,'%.2f');
        grid(ax,'on'); box(ax,'on');
        ylabel(ax,labels{coordinate},'Interpreter','latex','FontSize',22);
        low = min([entry.analytical(:,coordinate);entry.nonlinear(:,coordinate)]);
        high = max([entry.analytical(:,coordinate);entry.nonlinear(:,coordinate)]);
        pad = max((high-low)*0.18,1e-7);
        ylim(ax,[low-pad,high+pad]);
        if coordinate<3, ax.XTickLabel = []; end
    end
    xlabel(axes_list(3),'Time $t$ (-)','Interpreter','latex','FontSize',22);
    linkaxes(axes_list,'x');
    lgd = legend(axes_list(3),[h1,h2],{'Analytical cosine (linear mode)','Nonlinear RK4 simulation'}, ...
        'Interpreter','latex','FontSize',19,'Orientation','horizontal','Box','off');
    lgd.Layout.Tile = 'south';
    drawnow;
    stem = fullfile(output_dir,sprintf('mode_%d_comparison',mode));
    exportgraphics(fig,[stem,'.png'],'Resolution',300);
    savefig(fig,[stem,'.fig']);
end
