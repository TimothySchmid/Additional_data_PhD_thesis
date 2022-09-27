%
%--------------------------------------------------------------------------
% FUNCTION NAME:
%   plot_DVC_displacements
%
% DESCRIPTION:
%   Takes displacement fields from .mat structures obtained from DVC
%   analyses in DaVis. 2D plots can be made in three orthogonally
%   directions. 
%
% REQUIRED FUNCTIONS
%   - fct_colormap
%   - fct_plot_xslice
%   - fct_plot_yslice
%   - fct_plot_zslice
%   - fct_save_plots
%
% INPUT:
%
%   - INPUT - (structure) containing all input information:
%       - experimentname (str)  --> see name convention
%       - slice (str)           --> which direction (slice normal to given axis)
%       - slice_position (int)  --> position along axis
%       - displacement (str)    --> which displacement component
%           - 'x_slice'
%           - 'y_slice'
%           - 'z_slice'
%           - '3D_tot_displacement',
%           - '2D_tot_displacement' (tot displacement projected into plane)
%
%       - color_map_name (str)  --> for names: "help fct_colormap"
%       - contour_lines (str)   --> 'yes' or 'no' for activating/deactivating
%       - contour_color (str)   --> define color (rgb triplet, presets or hex)
%       - contour_width (float) --> line thickness
%       - vector_arrows (str)   --> 'yes' or 'no' for activating/deactivating
%       - vector_incr (int)     --> spacing of vector arrows w.r.t. grid
%       - vector_color (str)    --> define color (rgb triplet, presets or hex)
%       - vector_width (float)  --> line thickness
%       - print_png (str)       --> 'yes' or 'no' for activating/deactivating
%
%   - CELL  - (structure) containing all displacement fields from all time
%     steps, and all coordinates
%
% OUTPUT: 
%   - .png plots if needed
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

clear;
close all;
clc

% DEFINE INPUT
% ======================================================================= %

    INPUT.experimentname = 'EXP_823';
    INPUT.slice          = 'y_slice';
    INPUT.slice_position = 240;
    INPUT.displacement   = 'z_displacement';
    
    INPUT.color_map_name = 'roma';
    
    INPUT.contour_lines  = 'yes';
    INPUT.contour_color  = [0.1 0.1 0.1];
    INPUT.contour_width  = 1.2;
    
    INPUT.vector_arrows  = 'yes';
    INPUT.vector_incr    = 3;
    INPUT.vector_color   = 'red';
    INPUT.vector_width   = 1.0;
    
    INPUT.print_png      = 'yes';

% SET PATHS
% ======================================================================= %

    path_main =  pwd;
    path_data  = [pwd '/' INPUT.experimentname '_displacement_structure.mat'];
   
    path_colormap = [pwd '/_cmaps'];
    addpath(path_colormap)
    
% LOAD DISPLACEMENT STRUCTURE
% ======================================================================= %
    
    loadvar = path_data;
    load(loadvar)
         
% PLOTTING
% ======================================================================= %

switch INPUT.slice
    case 'x_slice'
        fct_plot_xslice(INPUT, CELL)
        
    case 'y_slice'
        fct_plot_yslice(INPUT, CELL)
        
    case 'z_slice'
        fct_plot_zslice(INPUT, CELL)
        
    otherwise
        error('slice name not defined')
end
