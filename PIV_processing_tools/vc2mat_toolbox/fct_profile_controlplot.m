function fct_profile_controlplot(cmap,X,Y,is_valid,TOPO,iRead,shift,INPUT,dt,plot_xcoords_profile,plot_ycoords_profile)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
valid_x = X(find(is_valid));
valid_y = Y(find(is_valid));

ch = convhull(valid_x,valid_y,'Simplify',true);

figure(1)
clf
set(gcf,'Units','normalized','Position',[.2 .2 .5 .6])
colormap(flipud(cmap))

s = fill(valid_x(ch),valid_y(ch),[0.5 0.5 0.5]);
hold on

pc = pcolor(X,Y,TOPO);
pc.FaceColor = 'interp';
pc.EdgeColor = 'none';

cb = colorbar;
ylabel(cb, 'Topography (mm)')
caxis([-5 5])

title(['Time: ' num2str((iRead+shift) * dt / 60) ' min'])
xlabel('Model width (mm)')
ylabel('Model length (mm)')

hAx=gca;
hAx.LineWidth=2;
hAx.FontSize = 14;

axis equal
axis([0 410 -470 0])
yticks(-450:50:0)
yticklabels({'450','400','350','300','250','200','150','100','50','0'})

box on
set(gca, 'Layer', 'Top')

plot(plot_xcoords_profile,plot_ycoords_profile,'r-.','LineWidth',2)
plot(valid_x(ch),valid_y(ch),'-k','LineWidth',3)

drawnow
end
