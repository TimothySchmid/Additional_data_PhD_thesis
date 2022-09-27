%
%--------------------------------------------------------------------------
% FUNCTION NAME:
%   load_DVC_displacements
%
% DESCRIPTION:
%   Takes displacement fields from .txt files (transformed from .dat files)
%   from te DVC analysis in DaVis and saves .mat structures containing all
%   timestep, displacement field components and coordinates
%
% INPUT:
%   - .txt files - Displacement fields from DaVis DVC analysis 
%
% OUTPUT: 
%   - CELL - (structure) Containing all displacement fields from all time
%     steps, and all coordinates
%
% ASSUMPTIONS AND LIMITATIONS:
%   None
% 
% FURTHER INFORMATION:
%
%  For more information, see <a href="matlab:
%  web('https://doi.org/10.1016/j.tecto.2021.229174')
%  ">Schmid et al., 2021</a>.
%
%  For more information, see <a href="matlab:
%  web('https://github.com/TimothySchmid/Characteristics_of_rotational_rifting.git')
%  ">Git hub repository</a>.
%
%  Latest DaVis readimx version for MacOS and Windows: <a href="matlab:
%  web('https://www.lavision.de/en/downloads/software/matlab_add_ons.php')
%  ">DaVis readimx</a>.
%--------------------------------------------------------------------------

% Author: Timothy Schmid, MSc., geology
% Institute of Geological Sciences, University of Bern
% Baltzerstrasse 1, Office 207
% 3012 Bern, CH
% email address: timothy.schmid@geo.unibe.ch
% January 2022; Last revision: 14/02/2022 
% Successfully tested on a Mac 64 bit using macOS Mojave
% (Vers. 10.14.6) and MATLABR2020b
%
%--------------------------------------------------------------------------

clear;
close all;
clc

% CHOOSE EXPERIMENT NAME INPUT
% ======================================================================= %
INPUT.experimentname = 'EXP_823';
INPUT.timestep       = 240;

% SET PATHS
% ======================================================================= %

    path_main =  pwd;
    path_data  = [pwd '/data_' INPUT.experimentname];
        
    path_colormap = [pwd '/_cmaps'];
    addpath(path_colormap)

% FIRST GET PROPER DIMENSIONS TO INITIALISE DATA STRUCTURE AND WRITE
% COORDINATES
% ======================================================================= %

    fprintf('\n make coordinate grid assemblage\n')
    cd(path_data)

    files = dir('*.txt');    l = length(files);

    name = files(1).name;
    disp_data = importdata(name);
    
    CELL.exp.name = INPUT.experimentname;
    CELL.exp.dt       = INPUT.timestep;

    % load coordinate vectors
    x = single(disp_data.data(:,1));
    y = single(disp_data.data(:,2));
    z = single(disp_data.data(:,3));
    
    % get proper dimensions by looking for unique values
    CELL.coords.xcoords = unique(x) - mean(unique(x));
    CELL.coords.ycoords = unique(y) - mean(unique(y));
    CELL.coords.zcoords = unique(z) - mean(unique(z));
    
    CELL.coords.xdim = uint16(numel(CELL.coords.xcoords));
    CELL.coords.ydim = uint16(numel(CELL.coords.ycoords));
    CELL.coords.zdim = uint16(numel(CELL.coords.zcoords));
    
    [CELL.coords.X,CELL.coords.Y,CELL.coords.Z] = ...
        ndgrid(CELL.coords.xcoords,CELL.coords.ycoords,CELL.coords.zcoords);
    
    clearvars disp_data
    
% GET DISPLACEMENTS
% ======================================================================= %
    
% load displacements
fprintf('\n make displacement assemblage\n')
    for istep = 1:l
        name = files(istep).name;
        disp_data = importdata(name);

        disp_x = single(disp_data.data(:,4));
        disp_y = single(disp_data.data(:,5));
        disp_z = single(disp_data.data(:,6));
        CELL.disp.X{istep} = fliplr(reshape(disp_x,CELL.coords.xdim,CELL.coords.ydim,CELL.coords.zdim));
        CELL.disp.Y{istep} = fliplr(reshape(disp_y,CELL.coords.xdim,CELL.coords.ydim,CELL.coords.zdim));
        CELL.disp.Z{istep} = fliplr(reshape(disp_z,CELL.coords.xdim,CELL.coords.ydim,CELL.coords.zdim));
        clearvars disp_data
    end
    
% calculate total displacement
fprintf('\n calculate total displacement\n')
    for istep = 1:l
        CELL.disp.TOT3D{istep} = sqrt(CELL.disp.X{istep}.^2 + ...
            CELL.disp.Y{istep}.^2 + CELL.disp.Z{istep}.^2);
    end
    
% SAVE DISPLACEMENTS
% ======================================================================= %
fprintf('\n save structure\n')
savevar = [path_data '/displacement_structure'];
save(savevar, 'CELL')
