function fct_vector_plot(nvec,COORDS,VEL)
    
% fct_vector_plot
%
% PLOT VECTORS
%
% Input:
%   - nvec     --> Vector spacing
%   - VEL      --> Velocity structure
%   - COORDS   --> Coodinates structure
%
% ======================================================================= %

        xvec    = round(size(COORDS.X,1)/nvec);
        yvec    = round(size(COORDS.X,2)/nvec);

        xv1     = round(linspace(2,size(COORDS.X,1)-1,xvec));
        xv2     = round(linspace(2,size(COORDS.Y,2)-1,yvec));

        
        q = quiver(COORDS.X(xv1,xv2),COORDS.Y(xv1,xv2),...
            VEL.U(xv1,xv2),VEL.V(xv1,xv2),...
            'color',[1.0, 0.2, 0.2],'linewidth',1.5);
        q.AutoScale = 'on';
        q.AutoScaleFactor = 0.5;
end

