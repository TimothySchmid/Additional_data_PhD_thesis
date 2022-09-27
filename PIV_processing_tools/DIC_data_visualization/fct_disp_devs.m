function [DISP,DISP_DEVS] = fct_disp_devs(DISP,COORDS,GRID)
        
% fct_disp_devs
%
% CALCULATES VELOCITY DERIVED VALUES
%
% Input:
%   - DISP      --> Displacement structure
%   - COORDS    --> Coodinates structure
%   - GRID      --> Structure containing grid information
%
% Output:
%   - DISP      --> Displacement structure
%   - DISP_DEVS --> Displacement derived values
%
% ======================================================================= %

% total displacements in 2d and 3d
    DISP.D2d = sqrt(DISP.Dx.^2 + DISP.Dy.^2);
    DISP.D3d = sqrt(DISP.Dx.^2 + DISP.Dy.^2 + DISP.Dz.^2);

 % Displacement gradients
    [Dxdy,Dxdx] = gradient(DISP.Dx,COORDS.dy*GRID.ny,COORDS.dx*GRID.nx); % calculate Dx/dy & Dx/dx
    [Dydy,Dydx] = gradient(DISP.Dy,COORDS.dy*GRID.ny,COORDS.dy*GRID.ny); % calculate Dy/dy & Dy/dx

 % Strain components
    DISP_DEVS.Exx = Dxdx;
    DISP_DEVS.Exy = Dxdy;
    DISP_DEVS.Eyx = Dydx;
    DISP_DEVS.Eyy = Dydy;
    
    DISP_DEVS.Gavg  =  (DISP_DEVS.Exy + DISP_DEVS.Eyx)/2;  % Avg. Shear
    DISP_DEVS.Div   = -(DISP_DEVS.Exx + DISP_DEVS.Eyy);    % Divergence in xy-plane
    
  % Invariants
    DISP_DEVS.I1 = DISP_DEVS.Exx  +  DISP_DEVS.Eyy;                                  % tr(S)
    DISP_DEVS.I2 = DISP_DEVS.Exx .*  DISP_DEVS.Eyy - DISP_DEVS.Exy .* DISP_DEVS.Eyx; % det(S) 
  
  % Angle of principal stretching axes 
    A = DISP_DEVS.Exx + DISP_DEVS.Eyy;
    B = DISP_DEVS.Exx - DISP_DEVS.Eyy;
    C = DISP_DEVS.Exy;
    
    DISP_DEVS.Theta_P = 0.5*atan(2*A ./ B);
    
  % Minimum principal stretching axes
    DISP_DEVS.Emin = 0.5*(A) - ((0.5*(B)).^2 + C.^2).^0.5;
    
  % Maximum principal stretching axes
    DISP_DEVS.Emax = 0.5*(A) + ((0.5*(B)).^2 + C.^2).^0.5;
    
  % max shear strain
    DISP_DEVS.Gmax = DISP_DEVS.Emax - DISP_DEVS.Emin;
end

