% plotting displacements and strains "on the fly"

clear
close all
clc


% INPUT
% ----------------------------------------------------------------------- %
  INPUT.experimentname      = 'R5';
  INPUT.datatype            = 'R5_ext_inc5_diff_11-3';
  
  INPUT.displacement_type   = 'cumulative';
  INPUT.strain_type         = 'infinitesimal';
% ----------------------------------------------------------------------- %


% SET PATHS TO DATA
% ----------------------------------------------------------------------- %
path_main = pwd;
path_to_experiment = [path_main, '/', INPUT.experimentname];
path_to_data = [path_to_experiment, '/', INPUT.datatype '_clean'];

addpath(path_main)
addpath 'customcolormap'
cmap = customcolormap([0 .25 .5 .75 1], ...
    {'#9d0142','#f66e45','#ffffbb','#65c0ae','#5e4f9f'});

% load metadata and assemble coordinate grid
cd([path_to_experiment, '/metadata'])
load('coordinate_system.mat')
load('meta_data.mat')

dt = EXP.dt;
% xcoords = xcoords - (-48.3577);
% ycoords = ycoords - 47.6794;

[X,Y]   = ndgrid(xcoords,ycoords);
dx      = diff(xcoords(1:2));
dy      = diff(ycoords(1:2));
[nx,ny] = size(X);

cd(path_to_data)
files = dir('*.mat');
n = length(files);

for iRead = n%(1:n)
    
file_now = files(iRead).name;

load(file_now)

% find square which contains no NaN's
  midx = round(size(is_valid,1)/2,0);
  midy = round(size(is_valid,2)/2,0);

% lines through mid point
  xline = is_valid(midx,:);
  yline = is_valid(:,midy);

% check for valid points along profiles
  xvalid = ~isnan(xline);
  yvalid = ~isnan(yline);
 
% get starting and end points of valid values
  ystart = find(xline,1,'first'); yend = find(xline,1,'last');
  xstart = find(yline,1,'first'); xend = find(yline,1,'last');
    
  xdir_data = xstart+10:xend-10;
  ydir_data = ystart+10:yend-10; 
  
% get clean data mask
  nx_clean = length(xdir_data);
  ny_clean = length(ydir_data);
  
  data_fill  = ones(nx_clean,ny_clean);
  whole_mask = NaN(size(X));
  
  whole_mask(xdir_data(1):xdir_data(end),ydir_data(1):ydir_data(end)) = data_fill;
  
  
% CALCULATE STRAINS
% ----------------------------------------------------------------------- %
% Displacement gradient tensor H
  [~,dudx] = gradient(Du, dx);
  [dudy,~] = gradient(Du, dy);
  [~,dvdx] = gradient(Dv, dx);
  [dvdy,~] = gradient(Dv, dy);

% Deformation gradient tensor F
  F(1,1,:,:) = dudx + 1;
  F(1,2,:,:) = dudy;
  F(2,1,:,:) = dvdx;
  F(2,2,:,:) = dvdy + 1;

% Calculate strains
  PLT = fct_calculate_strain(INPUT,nx,ny,F,whole_mask);

plot_what = PLT.eyy;

figure(1)
clf
set(gcf,'Units','normalized','Position',[.2 .2 .5 .6])
colormap(cmap)

thresh_map = (plot_what.*whole_mask);
thresh_map(thresh_map >= -0.01) = 1;

switch INPUT.displacement_type
    case 'incremental'
%         pcolor(X,Y,plot_what.*whole_mask)
        pcolor(X,Y,plot_what)
    case 'cumulative'
        pcolor(X+Du,Y-Dv,plot_what.*whole_mask)
    otherwise
end

shading flat

axis equal
axis([-250 250 -300 250])

colorbar
% caxis([-0.1 0.1])
drawnow


end

restoredefaultpath
cd(path_main)



