function figure_format(fig_handle)
% Create axes, modify the text size and plot size
axes1 = axes('Parent',fig_handle,'FontSize',14,...
    'Position',[0.0551 0.0809 0.9185 0.873]);
box(axes1,'on');
hold(axes1,'on');
set(gcf,'Position',[100 100 1140 600])

