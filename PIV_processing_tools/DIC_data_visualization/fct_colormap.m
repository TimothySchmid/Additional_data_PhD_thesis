function [c_map] = fct_colormap(INPUT)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
switch INPUT.color_map_name
    
    case 'viridis'
        c = load('viridis.mat');
        c_map = c.cm;
        
    case 'inferno'
        c = load('inferno.mat');
        c_map = c.cm;
        
    case 'plasma'
        c = load('plasma.mat');
        c_map = c.cm;
        
    case 'magma'
        c = load('magma.mat');
        c_map = c.cm;
        
    case 'vik'
        c = load('vik.mat');
        c_map = c(:,:).vik;
        
    case 'roma'
        c = load('roma.mat');
        c_map = c(:,:).roma;
        
    case 'davos'
        c = load('davos.mat');
        c_map = c(:,:).davos;
        
    case 'oleron'
        c = load('oleron.mat');
        c_map = c(:,:).oleron;
        
    case 'cork'
        c = load('cork.mat');
        c_map = c(:,:).cork;
        
    case 'broc'
        c = load('broc.mat');
        c_map = c(:,:).broc;
        
    case 'vik10'
        c = load('vik10.mat');
        c_map = c(:,:).vik10;
        
    case 'vik25'
        c = load('vik25.mat');
        c_map = c(:,:).vik25;
        
    case 'rgb'
        
        c_map = [ ...
            94    79   162
            50   136   189
            102   194   165
            171   221   164
            230   245   152
            255   255   191
            254   224   139
            253   174    97
            244   109    67
            213    62    79
            158     1    66  ] / 255;
end
    
end

