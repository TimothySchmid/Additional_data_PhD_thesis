% plotting topographic profiles and displacements

clear
close all
clc


% INPUT
% ----------------------------------------------------------------------- %
  INPUT.experimentname      = 'R1';
  INPUT.datatype            = 'R1_ext_inc6_diff_15-5';
  INPUT.displacement        = 'R1_ext_inc6_sumdiff_15-5';
  
  INPUT.profile_pos         = [100, 200, 300];
  INPUT.control_plot        = 'no';
  
  INPUT.colormap            = 'roma';
  INPUT.save                = 'yes';
% ----------------------------------------------------------------------- %


% SET PATHS TO DATA
% ----------------------------------------------------------------------- %
  path_main = pwd;
  path_to_experiment = [path_main, '/', INPUT.experimentname];
  path_to_data = [path_to_experiment, '/', INPUT.datatype '_clean'];
  path_to_png = [path_to_experiment, '/png_profiles'];
  path_to_disp = [path_to_experiment, '/', INPUT.displacement '_clean'];
  
  addpath(path_main)

  addpath([path_main, '/_cmaps'])
  cmap = fct_colormap(INPUT);

  cd(path_to_experiment)
  mkdir png_profiles
  
% load metadata and assemble coordinate grid
  cd([path_to_data, '/metadata'])
  load('coordinate_system.mat')
  load('meta_data.mat')
  dt = 3*600;%EXP.dt;
  
  xcoords = xcoords - mean(xcoords);
  ycoords = ycoords - mean(ycoords);
  [X,Y]   = ndgrid(xcoords,ycoords);
  [nx,ny] = size(X);

  cd(path_to_data)
  files = dir('*.mat');
  n = length(files);
  
% ASSEMBLE STRUCTURE FOR CSV FILE
% ----------------------------------------------------------------------- %
profile_1 = struct('Position',{},...
                   'Coordinates',{},...
                   'Topography',{});
               
profile_2 = struct('Position',{},...
                   'Coordinates',{},...
                   'Topography',{});
               
profile_3 = struct('Position',{},...
                   'Coordinates',{},...
                   'Topography',{});

DATA = struct('Time',{},'Profile_1',profile_1,'Profile_2',profile_2,'Profile_3',profile_3);

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

clearvars Du Dv Dw

cd(path_to_disp)
files_disp = dir('*.mat');
disp_now = files_disp(iRead+1).name;
load(disp_now,'Du','Dv','Dw');
cd(path_to_data)

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

% PROFILES
% ----------------------------------------------------------------------- %
profile_num = 1;
for profile_pos = INPUT.profile_pos
% find positions for profile
[res, profile_ind] = min(abs(xcoords-profile_pos));

TOPO_profile_all = TOPO(profile_ind,:);
profile_mask = ~isnan(TOPO_profile_all);

topo_pos_first = find(profile_mask, 1, 'first');
topo_pos_last  = find(profile_mask, 1, 'last');

plot_ycoords_profile = (ycoords(topo_pos_first:topo_pos_last));
plot_topo_profile    = (TOPO(profile_ind,topo_pos_first:topo_pos_last));

plot_xcoords_profile = xcoords(profile_ind) * ones(size(plot_ycoords_profile));

switch INPUT.control_plot
    case 'yes'
        fct_profile_controlplot(cmap,X,Y,is_valid,TOPO,iRead,shift,INPUT,...
            dt,plot_xcoords_profile,plot_ycoords_profile)
    case 'no'
        figure(2)
        clf
        set(gcf,'Units','normalized','Position',[.1 .3 .8 .3])
        
        % topography
        yl = yline(0,'-.','LineWidth',2);
        yl.Color = [0.5 0.5 0.5];
        hold on
        plot(plot_ycoords_profile,plot_topo_profile,'LineWidth',3)
        
        %lateral boundaries
        xleft  = plot_ycoords_profile(1);
        xright = plot_ycoords_profile(end);
        ytop   = 5;
        ybot   = -ytop;
        
        plot([xleft xleft],[ybot, ytop],'r','LineWidth',3)
        plot([xright xright],[ybot, ytop],'r','LineWidth',3)
        
        set (gca, 'xdir', 'reverse' )
        daspect([3, 1, 1])
        axis([-470, 0, -10, 10])
        
        xticks(-450:50:0)
        xticklabels({'450','400','350','300','250','200','150','100','50','0'})
        xlabel('Model width (mm)','FontSize', 16)
        ylabel('Topography (mm)')
        
        hAx=gca;
        hAx.LineWidth=2;
        hAx.FontSize = 14;
        
        box on
        set(gca, 'Layer', 'Top')
             
        title(['Time: ' num2str((iRead+shift) * dt / 60) ' min'],'FontSize',18)
        
    otherwise
        error('Unclear if control plot is requested. Check spelling')
end
        
% SAVING PLOTS AS PNG
% ----------------------------------------------------------------------- %
switch INPUT.save
    case 'yes'
        exportgraphics(hAx,[path_to_png,'/',INPUT.datatype,'_pos_',...
            num2str(profile_pos),'_',num2str((iRead+shift) * dt / 60),...
            '_min.png'],'Resolution',300)
        
        saveas(gcf,[path_to_png,'/',INPUT.datatype,'_pos_',...
            num2str(profile_pos),'_',num2str((iRead+shift) * dt / 60),...
            '_min'],'epsc')
    case 'no'
    otherwise
        error('unclear if saving is requested. Check spelling')
end

% WRITE DATA TO STRUCTURE
% ----------------------------------------------------------------------- %
DATA{iRead}.Time = (iRead+shift) * dt / 60;

if profile_num == 1
    DATA{iRead}.profile_1.Position    = INPUT.profile_pos(1);
    DATA{iRead}.profile_1.Coordinates = plot_ycoords_profile;
    DATA{iRead}.profile_1.Topography  = plot_topo_profile;
elseif profile_num == 2
    DATA{iRead}.profile_2.Position    = INPUT.profile_pos(2);
    DATA{iRead}.profile_2.Coordinates = plot_ycoords_profile;
    DATA{iRead}.profile_2.Topography  = plot_topo_profile;
else
    DATA{iRead}.profile_3.Position    = INPUT.profile_pos(3);
    DATA{iRead}.profile_3.Coordinates = plot_ycoords_profile;
    DATA{iRead}.profile_3.Topography  = plot_topo_profile;
end

profile_num = profile_num + 1;

end    
end

% SAVE STRUCTURE
% ----------------------------------------------------------------------- %
cd(path_to_experiment)
if contains (INPUT.datatype,'short')
    savevar = 'profile_data_shortening';
    save(savevar,'DATA');
else
    savevar = 'profile_data_extension';
    save(savevar,'DATA');
end

restoredefaultpath
cd(path_main)
