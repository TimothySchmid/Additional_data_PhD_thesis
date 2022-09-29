% plotting displacements and strains "on the fly"

clear
close all
clc


% INPUT
% ----------------------------------------------------------------------- %
  INPUT.experimentname      = 'R5';
  INPUT.datatype            = 'R5_ext_inc5_diff_11-3';
  INPUT.plot_style          = 'topview';
  INPUT.colormap            = 'roma';
  INPUT.save                = 'no';
% ----------------------------------------------------------------------- %


% SET PATHS TO DATA
% ----------------------------------------------------------------------- %
  path_main = pwd;
  path_to_experiment = [path_main, '/', INPUT.experimentname];
  path_to_data = [path_to_experiment, '/', INPUT.datatype '_clean'];
  path_to_png = [path_to_experiment, '/png_topography'];
  
  addpath(path_main)

  addpath([path_main, '/_cmaps'])
  cmap = fct_colormap(INPUT);

  cd(path_to_experiment)
  mkdir png_topography
  
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
      shift = 0; %n;
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

TOPO(1:xind_first + cut ,:) = NaN;
TOPO(xind_last - cut:end,:) = NaN;
TOPO(:,1:yind_first + cut)  = NaN;
TOPO(:,yind_last - cut:end) = NaN;

valid_x = X(find(is_valid));
valid_y = Y(find(is_valid));

ch = convhull(valid_x,valid_y,'Simplify',true);

switch INPUT.plot_style
    case 'topview'
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
        caxis([-1.5 1.5])
        
        title(['Time: ' num2str((iRead+shift) * dt / 60) ' min'])
        xlabel('Model width (mm)')
        ylabel('Model length (mm)')
        
        hAx=gca;
        hAx.LineWidth=2;
        hAx.FontSize = 14;
        
        axis equal
        axis([0 400 -450 0])
        yticks(-450:50:0)
        yticklabels({'450','400','350','300','250','200','150','100','50','0'})
        
        box on
        set(gca, 'Layer', 'Top')
        
        plot(valid_x(ch),valid_y(ch),'-k','LineWidth',3)
        
%         yline(0,'-r','LineWidth',3)
%         xline(0,'-r','LineWidth',3)
        drawnow
        
    case 'shading'
        figure(2)
        clf
        set(gcf,'Units','normalized','Position',[.2 .2 .5 .6])
        colormap(flipud(cmap))
        
        z_comp = -4 * ones(length(valid_x(ch)));
        s = fill3(valid_x(ch),valid_y(ch),z_comp,[0.5 0.5 0.5]);
        hold on
        
        sf = surf(X,Y,TOPO);
        sf.FaceColor = 'interp';
        sf.EdgeColor = 'none';
        hold on
        
        grid off
        view(2)
        
        lightangle(45,85)
        sf.FaceLighting = 'gouraud';
        sf.AmbientStrength = 0.8;
        sf.DiffuseStrength = 0.8;
        sf.SpecularStrength = 0.9;
        sf.SpecularExponent = 25;
        sf.BackFaceLighting = 'unlit';
        
        cb = colorbar;
        cb.FontSize = 16;
        ylabel(cb, 'Topography (mm)')
        caxis([-4.5 4.5])
        
        title(['Time: ' num2str((iRead+shift) * dt / 60) ' min'])
        xlabel('Model width (mm)')
        ylabel('Model lengtht (mm)')
        
        hAx=gca;
        hAx.LineWidth=2;
        hAx.FontSize = 16;
        
        axis equal
        axis([0 400 -450 0])
        yticks(-450:50:0)
        yticklabels({'450','400','350','300','250','200','150','100','50','0'})
        
        box on
        set(gca, 'Layer', 'Top')
        
        height_val = 3*ones(length(valid_x(ch)));
        plot3(valid_x(ch),valid_y(ch),height_val,'-k','LineWidth',3)
        drawnow
        
    otherwise
        error('unclear which plot option is requested. Check spelling')
end

% SAVING PLOTS AS PNG
% ----------------------------------------------------------------------- %
switch INPUT.save
    case 'yes'
        exportgraphics(hAx,[path_to_png,'/',INPUT.datatype,'_topogrpahy_',...
            num2str((iRead+shift) * dt / 60),' min.png'],'Resolution',600)
    case 'no'
    otherwise
        error('unclear if saving is requested. Check spelling')
end
    
end

restoredefaultpath
cd(path_main)
