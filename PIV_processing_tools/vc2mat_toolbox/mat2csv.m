% create .csv files

clear
close all
clc


% INPUT
% ----------------------------------------------------------------------- %
  INPUT.experimentname      = 'R5';
% ----------------------------------------------------------------------- %


% SET PATHS TO DATA
% ----------------------------------------------------------------------- %
  path_main = pwd;
  path_to_experiment = [path_main, '/', INPUT.experimentname];

  cd(path_to_experiment)
  
  
% WRITE EXTENSION PHASE .CSV FILE
% ----------------------------------------------------------------------- %
  loadvar = 'profile_data_extension.mat';
  load(loadvar)

  n = size(DATA,2);
  header_cell = cell(1,2*n);
  
  % go through data, get data per time step and save it in cell
  pos = 1;
  for iRead = (1:n)
      time_now        =  DATA{iRead}.Time;
      position_now    =  DATA{iRead}.profile_3.Position;
      coordinates_now = -DATA{iRead}.profile_3.Coordinates;
      topography_now  = (DATA{iRead}.profile_3.Topography);
      
      TIME{iRead} = time_now;
      CELL{pos}   = coordinates_now';
      CELL{pos+1} = topography_now';
      
      header_cell(pos)   = {[num2str(time_now),' min Coordinates [mm]']};
      header_cell(pos+1) = {[num2str(time_now),' min Topography [mm]']};
      
      pos = pos + 2;
  end

  % find maximum column size for homogenisation
    max_array = zeros(1,2*n);
    for cellpos = (1:2*n)
       max_array(1,cellpos) = size(CELL{cellpos},1); 
    end

  % make homogenized matrix filled with NaN values and fill in columns
    table_matrix = NaN(max(max_array),2*n);
    for cellpos = (1:2*n)
        col_now = CELL{cellpos};
       table_matrix(1:length(col_now),cellpos) = col_now; 
    end

  % write table and make .csv file
    header_cell = string(header_cell);
    T = array2table(table_matrix,'VariableNames',header_cell);
    writetable(T,[INPUT.experimentname, '_extension_profile_3.csv'])


% WRITE SHORTENING PHASE .CSV FILE
% ----------------------------------------------------------------------- %
  loadvar = 'profile_data_shortening.mat';
  load(loadvar)

  n = size(DATA,2);
  header_cell = cell(1,2*n);
  
  % go through data, get data per time step and save it in cell
  pos = 1;
  for iRead = (1:n)
      time_now        =  DATA{iRead}.Time;
      position_now    =  DATA{iRead}.profile_3.Position;
      coordinates_now = -DATA{iRead}.profile_3.Coordinates;
      topography_now  = (DATA{iRead}.profile_3.Topography);
      
      TIME{iRead} = time_now;
      CELL{pos}   = coordinates_now';
      CELL{pos+1} = topography_now';
      
      header_cell(pos)   = {[num2str(time_now),' min Coordinates [mm]']};
      header_cell(pos+1) = {[num2str(time_now),' min Topography [mm]']};
      
      pos = pos + 2;
  end

  % find maximum column size for homogenisation
    max_array = zeros(1,2*n);
    for cellpos = (1:2*n)
       max_array(1,cellpos) = size(CELL{cellpos},1); 
    end

  % make homogenized matrix filled with NaN values and fill in columns
    table_matrix = NaN(max(max_array),2*n);
    for cellpos = (1:2*n)
        col_now = CELL{cellpos};
       table_matrix(1:length(col_now),cellpos) = col_now; 
    end

  % write table and make .csv file
    header_cell = string(header_cell);
    T = array2table(table_matrix,'VariableNames',header_cell);
    writetable(T,[INPUT.experimentname, '_shortening_profile_3.csv'])
