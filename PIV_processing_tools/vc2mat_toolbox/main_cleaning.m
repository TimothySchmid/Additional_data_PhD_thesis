clear
close all
clc


% INPUT PARAMETERS
% ======================================================================= %
EXP.name              = 'R5';                         % --> 'experiment name'
EXP.data_choice       = 'R5_ext_inc5_diff_11-3';   % --> 'data set name'


EXP.check_plot        = 'yes';                         % --> 'yes'/'no'
EXP.cleaning_type     = 'DCT-PLS';                    % --> 'default'/'DCT-PLS'

% for default cleaning
EXP.outlier.threshmed = 1;          %1                % --> noise threshold
EXP.outlier.eps       = 1e-1;       %1e-2             % --> estimated measurement noise level
EXP.outlier.neighbour = 7;                            % --> neighborhood radius: 1 = 3x3, 2 = 5x5 etc.

% for DCT-PLS algorithm
EXP.DCT_PLS           = 1e-2;                         % --> value for DCT-PLS data correction
% ======================================================================= %


% SET PATHS TO FUNCTIONS AND MAKE FOLDERS
% ----------------------------------------------------------------------- %
parent_path  = pwd;
addpath(parent_path)
addpath 'customcolormap'

if isunix
    addpath([parent_path, '/readimx-v2.1.8-osx']); % Mac
else
    addpath([parent_path, '/readimx-v2.1.9']);     % Windows
end

path_to_experiment = [parent_path, '/', EXP.name];
path_to_data = [path_to_experiment, '/', EXP.data_choice];

cd(path_to_experiment)

mkdir([EXP.data_choice, '_clean'])
clean_data_path = [path_to_experiment '/' EXP.data_choice '_clean'];

cd(clean_data_path)
mkdir metadata
metadata_path = [clean_data_path '/metadata'];

% get data list and make new folder for mat files
cd(path_to_data)
files = dir('*.vc7');
files(strncmp({files.name}, '.', 1)) = [];
n = length(files);


% LOCATE DISPLACEMENT COMPONENTS
% ----------------------------------------------------------------------- %
vc_struc_init = readimx(files(1).name);

% search for correct places
loc_u = fct_find_location(vc_struc_init,'U0');
loc_v = fct_find_location(vc_struc_init,'V0');
loc_w = fct_find_location(vc_struc_init,'W0');
loc_h = fct_find_location(vc_struc_init,'TS:Height');
loc_m = fct_find_location(vc_struc_init,'MASK');


% SCALING VALUES FOR COORDINATE SYSTEM
% ----------------------------------------------------------------------- %
[slope_x, offset_x, step_x] = fct_get_scaling(vc_struc_init, 'X');
[slope_y, offset_y, step_y] = fct_get_scaling(vc_struc_init, 'Y');
[slope_z, offset_z, step_z] = fct_get_scaling(vc_struc_init, 'Z');
[slope_i, offset_i] = fct_get_scaling(vc_struc_init, 'I');


% ASSEMBLE COORDINATE SYSTEM
% ----------------------------------------------------------------------- %
dim = size(vc_struc_init.Frames{1}.Components{loc_u}.Planes{:});

xcoords = slope_x * (linspace(0, dim(2), dim(1)) * step_x) + offset_x;
ycoords = slope_y * (linspace(0, dim(1), dim(2)) * step_y) + offset_y;

% write experiment data and coordinates
savevar = [metadata_path '/coordinate_system'];
save(savevar, '*coords', 'slope*', 'offset*', 'step*', '-v7.3')


% WRITE META DATA
% ----------------------------------------------------------------------- %
EXP = fct_write_metadata(vc_struc_init, EXP);

if contains(EXP.data_choice,'sum')
   EXP.dt = 600; 
end

% clearvars vc_struc_init slope_* offset_* *coords dim savevar step_* ...
%     -except slope_i


% RUN THROUGH FILES
% ----------------------------------------------------------------------- %
tStart = tic;

fct_print_statement('starting')
fct_print_statement('cleaning')

for iRead = progress(1:n)
    
  % get step
    cd(path_to_data)
    step_now = files(iRead).name;

  % get current .vc7 structure
    vc_struc = readimx(step_now);
    
  % prepare mask
    M0 = logical(vc_struc.Frames{1}.Components{loc_m}.Planes{:});
    is_valid = fct_prepare_mask(M0);
    
    clean_mask = double(is_valid);
    clean_mask(clean_mask==0) = NaN;
    
  % get needed components
    U0 = vc_struc.Frames{1}.Components{loc_u}.Planes{:};
    V0 = vc_struc.Frames{1}.Components{loc_v}.Planes{:};
    W0 = vc_struc.Frames{1}.Components{loc_w}.Planes{:};
    H0 = vc_struc.Frames{1}.Components{loc_h}.Planes{:};
    
  % clean extracted buffers
    switch EXP.cleaning_type
        case 'default'
            % default (slower but more robust for "nasty" data)
            [H_temp, U_temp, V_temp, W_temp] = fct_clean_raw_data(H0, ...
                U0, V0, W0, EXP, is_valid);
        case 'DCT-PLS'
            % DCT-PLS filtering (works well and fast for "good" data)
            U_temp = fct_clean_raw_data_DCT(U0, clean_mask, is_valid, EXP);
            V_temp = fct_clean_raw_data_DCT(V0, clean_mask, is_valid, EXP);
            W_temp = fct_clean_raw_data_DCT(W0, clean_mask, is_valid, EXP);
            H_temp = fct_clean_raw_data_DCT(H0, clean_mask, is_valid, EXP);
        otherwise
            error('not clear how data should be cleared. Check spelling of EXP.cleaning_type')  
    end
    
  % height correction (if incremental data)
    if iRead == 1
        [correction_plane, boundaries] = fct_extract_data(H_temp, is_valid);
        [Dev, fit_vals]                = fct_correct_height(correction_plane);
         Dev_ext                       = fct_reassign_values(Dev, boundaries, H_temp);
         EXP.HeightCoefficients        = fit_vals;
         savevar = [metadata_path '/meta_data'];
         save(savevar, 'EXP', '-v7.3')
    end
    
  % scale new variables
    Du = single(U_temp * slope_i);
    Dv = single(V_temp * slope_i);
    Dw = single(W_temp * slope_i);
    TOPO = single(H_temp - Dev_ext) * slope_i;
    
  % control plot
    fct_check_plot(EXP, H0, U0, V0, W0, TOPO, Du, Dv, Dw, iRead)
    
  % write new data as 7.3 .mat file (can also be read in python)
    file_name = strrep(step_now, 'vc7', 'mat');
    savevar = [clean_data_path '/' file_name];
    
    if contains(EXP.data_choice,'sum')
        save(savevar, 'Du', 'Dv', 'Dw', 'is_valid', '-v7.3')
    else
        save(savevar, 'Du', 'Dv', 'Dw', 'TOPO', 'is_valid','-v7.3')
    end
    
end

fct_print_statement('ending', tStart)

restoredefaultpath
