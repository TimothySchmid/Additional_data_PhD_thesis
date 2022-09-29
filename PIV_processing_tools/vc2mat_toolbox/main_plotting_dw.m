% plotting displacements and strains "on the fly"

clear
close all
clc


% INPUT
% ----------------------------------------------------------------------- %
  INPUT.experimentname      = 'R1';
  INPUT.datatype            = 'R1_ext_inc6_diff_15-5';
  INPUT.colormap            = 'vik';
  INPUT.save                = 'yes'; 
  INPUT.xvectors            = 10;
% ----------------------------------------------------------------------- %


% SET PATHS TO DATA
% ----------------------------------------------------------------------- %
  path_main = pwd;
  path_to_experiment = [path_main, '/', INPUT.experimentname];
  path_to_data = [path_to_experiment, '/', INPUT.datatype '_clean'];
  path_to_png = [path_to_experiment, '/png_dw'];
  
  addpath(path_main)
  
  addpath([path_main, '/_cmaps'])
  cmap = fct_colormap(INPUT);

  cd(path_to_experiment)
  mkdir png_dw
  
% load metadata and assemble coordinate grid
  cd([path_to_data, '/metadata'])
  load('coordinate_system.mat')
  load('meta_data.mat')
  dt = 1800;%EXP.dt;
  
  xcoords = xcoords - mean(xcoords);
  ycoords = ycoords - mean(ycoords);
  [X,Y]   = ndgrid(xcoords,ycoords);
  [nx,ny] = size(X);

  cd(path_to_data)
  files = dir('*.mat');
  n = length(files);

% GET TIME SHIFT FOR SHORTENING PHASE
% ----------------------------------------------------------------------- %
  if contains (INPUT.datatype,'short')
      shift = 0;%n;
  else
      shift = 0;
  end
  
for iRead = (1:n)
    
file_now = files(iRead).name;
load(file_now)
 
if iRead == 1
   % check for proper position of data   
   y_pos_slice = round(nx/2);
   x_pos_slice = round(ny/2);
   
   y_array = is_valid(y_pos_slice,:);
   x_array = is_valid(:,x_pos_slice)';
   
   xind_first = find(x_array, 1, 'first');
   xind_last  = find(x_array, 1, 'last');
    
   yind_first = find(y_array, 1, 'first');
   yind_last  = find(y_array, 1, 'last');
      
   xcoords = xcoords - xcoords(xind_first);
   ycoords = ycoords - ycoords(yind_first);
   
   [X,Y] = ndgrid(xcoords,ycoords);
else
    y_pos_slice = round(nx/2);
    x_pos_slice = round(ny/2);
    
    y_array = is_valid(y_pos_slice,:);
    x_array = is_valid(:,x_pos_slice)';
    
    xind_first = find(x_array, 1, 'first');
    xind_last  = find(x_array, 1, 'last');
    
    yind_first = find(y_array, 1, 'first');
    yind_last  = find(y_array, 1, 'last');
end

% cut nice data
cut = 15;

is_valid(1:xind_first + cut ,:) = NaN;
is_valid(xind_last - cut:end,:) = NaN;
is_valid(:,1:yind_first + cut)  = NaN;
is_valid(:,yind_last - cut:end) = NaN;

Du(1:xind_first + cut ,:) = NaN;
Du(xind_last - cut:end,:) = NaN;
Du(:,1:yind_first + cut)  = NaN;
Du(:,yind_last - cut:end) = NaN;

Dv(1:xind_first + cut ,:) = NaN;
Dv(xind_last - cut:end,:) = NaN;
Dv(:,1:yind_first + cut)  = NaN;
Dv(:,yind_last - cut:end) = NaN;

Dw(1:xind_first + cut ,:) = NaN;
Dw(xind_last - cut:end,:) = NaN;
Dw(:,1:yind_first + cut)  = NaN;
Dw(:,yind_last - cut:end) = NaN;

valid_x = X(find(is_valid));
valid_y = Y(find(is_valid));

valid_x(isnan(valid_x)) = [];
valid_y(isnan(valid_y)) = [];

ch = convhull(valid_x,valid_y,'Simplify',true);

    figure(1)
    clf
    set(gcf,'Units','normalized','Position',[.2 .2 .5 .6])
    colormap(cmap)
    
    s = fill(valid_x(ch),valid_y(ch),[0.5 0.5 0.5]);
    hold on
    
    pc = pcolor(X,Y,Dw);
    pc.FaceColor = 'interp';
    pc.EdgeColor = 'none';
    
    cb = colorbar;
    cb.FontSize = 16;
    ylabel(cb, 'Displacement dw (mm)')
    caxis([-0.5 0.5])
    
    title(['Time: ' num2str((iRead+shift) * dt / 60) ' min'])
    xlabel('Model width (mm)')
    ylabel('Model length (mm)')
    
    hAx=gca;
    hAx.LineWidth=2;
    hAx.FontSize = 14;
    
    axis equal
    axis([0 400 -470 0])
    yticks(-450:50:0)
    yticklabels({'450','400','350','300','250','200','150','100','50','0'})
    
    box on
    set(gca, 'Layer', 'Top')
    
    plot(valid_x(ch),valid_y(ch),'-k','LineWidth',3)
    
% GET NICE VECTOR DISTRIBUTION
% ----------------------------------------------------------------------- %
% get length ratio of data
[x_size,y_size] = size(Du); mid_x = round(x_size/2); mid_y = round(y_size/2);
x_data_width = Du(:,mid_x); y_data_width = Du(mid_y,:);

x_valid_pos = find(~isnan(x_data_width)); lx = length(x_valid_pos);
y_valid_pos = find(~isnan(y_data_width)); ly = length(y_valid_pos);

x_vec_first = x_valid_pos(1); x_vec_last = x_valid_pos(end);
y_vec_first = y_valid_pos(1); y_vec_last = y_valid_pos(end);

ratio_yx_vec = ly/lx;   INPUT.yvectors = round(INPUT.xvectors * ratio_yx_vec);

margin = 10;

x_vec_array = round(linspace(x_vec_first+margin,x_vec_last-margin,INPUT.xvectors));
y_vec_array = round(linspace(y_vec_first+margin,y_vec_last-margin,INPUT.yvectors));

q = quiver(X(x_vec_array,y_vec_array),Y(x_vec_array,y_vec_array),...
           Du(x_vec_array,y_vec_array),-Dv(x_vec_array,y_vec_array),...
           'color',[0.2 0.2 0.2],'LineWidth',2);
q.AutoScaleFactor = 0.5;
    
%     yline(0,'-r','LineWidth',3)
%     xline(0,'-r','LineWidth',3)
    drawnow
    
% SAVING PLOTS AS PNG
% ----------------------------------------------------------------------- %
switch INPUT.save
    case 'yes'
        exportgraphics(hAx,[path_to_png,'/',INPUT.datatype,'_dw_',...
            num2str((iRead+shift) * dt / 60),' min.png'],'Resolution',600)
    case 'no'
    otherwise
        error('unclear if saving is requested. Check spelling')
end

end

restoredefaultpath
cd(path_main)
