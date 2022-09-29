function fct_print_statement(statement, varargin)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

switch statement
    case 'starting'
        fprintf('\n\nStart processing data -------------------------\n')
        fprintf('-----------------------------------------------\n')
        fprintf('\n \n')
        
    case 'cleaning'
        fprintf('\n\nClean and prepare data ------------------------\n')
        fprintf('-----------------------------------------------\n')
        fprintf('\n \n')
        
    case 'refinement'
        fprintf('\n\nInterpolate displacement on refined grid ------\n')
        fprintf('-----------------------------------------------\n')
        fprintf('\n \n')
        
    case 'material'
        fprintf('\n\nCalculate Lagrangian sum ----------------------\n')
        fprintf('-----------------------------------------------\n')
        fprintf('\n \n')
        
    case 'H'
        fprintf('\n\nCalculate displacement gradient tensor H ------\n')
        fprintf('-----------------------------------------------\n')
        fprintf('\n \n')
        
    case 'F'
        fprintf('\n\nCalculate deformation gradient tensor F -------\n')
        fprintf('-----------------------------------------------\n')
        fprintf('\n \n')
        
    case 'save'
        fprintf('\n\nSaving ----------------------------------------\n')
        fprintf('-----------------------------------------------\n')
        fprintf('\n \n')
        
    case 'summation'
        fprintf('\n\nMaterial displacement summation ---------------\n')
        fprintf('-----------------------------------------------\n')
        fprintf('\n \n')
        
    case 'ending'
        if nargin > 2
            varargin = tStart;
        end
        
        tEnd = toc(cell2mat(varargin));
        d = datetime;
        s1 = char(d);

        fprintf('\n')
        fprintf('-----------------------------------------------\n')
        fprintf(['Wall clock time: ', num2str(round(tEnd/60,2)), ' minutes\n'])
        fprintf('-----------------------------------------------\n')
        fprintf(['Processing finished on: ', s1, '\n'])
        fprintf('-----------------------------------------------\n')
        
    otherwise
end

end

