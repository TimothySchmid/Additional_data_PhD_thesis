function     DataC = fct_clean_raw_data_DCT(Data0, clean_mask, is_valid, EXP)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

dct_val = EXP.DCT_PLS;

% apply clean mask to crop data
  Data = Data0 .* clean_mask;
  
% turn remaining zero values (inside) to NaN values
  Data(Data==0) = NaN;
  
% fill surrounding space with zeros (must not be NaN)
  Data(is_valid==0) = 0;

% apply DCT-PLS filtering
  Data = smoothn(Data,dct_val);
  
% turn surrounding values back to NaN
  Data(is_valid==0) = NaN;
  
% return clean data
  DataC = Data;
end

