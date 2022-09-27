function fct_save_plots(INPUT, CELL, step)

% fct_save_plots
%
% PLOTTING SLICES NORMAL TO Z AXIS
%
% Input:
%   - INPUT      --> Structure containing initial input parameters
%   - CELL       --> Structure containing  displacements and coordinates
%   - step       --> file number to calculate absolute time
%
% ======================================================================= %

    switch INPUT.print_png
        case 'yes'
            if isfolder('png_figures')
                % pass
            else
                mkdir 'png_figures'
            end
            
            
            fig_name = [INPUT.experimentname, '_', INPUT.slice, '_pos_', ...
                num2str(INPUT.slice_position), '_', INPUT.displacement, '_time_',...
                num2str(step * CELL.exp.dt,'%4.4d')];
            
            print('-dpng','-r600','-noui',[pwd '/png_figures/' fig_name '.png'])
        case 'no'
            
        otherwise
            error('Unclear if figure should be printed or not')
    end
end
