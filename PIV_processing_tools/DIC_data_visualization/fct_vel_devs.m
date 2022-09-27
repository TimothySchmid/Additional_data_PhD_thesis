function [VEL,VEL_DEVS] = fct_vel_devs(DISP,COORDS,GRID,INPUT)

% fct_vel_devs
%
% CALCULATES VELOCITY DERIVED VALUES
%
% Input:
%   - DISP     --> Displacement structure
%   - COORDS   --> Coodinates structure
%   - GRID     --> Structure containing grid information
%   - INPUT    --> Structure containing initial input parameters
%
% Output:
%   - VEL      --> Velocity components
%   - VEL_DEVS --> Velocity derivative values
%
% ======================================================================= %

switch INPUT.disp_type
    case 'incremental'
        time_increment = INPUT.image_interval;
    case 'cumulative'
        time_increment = GRID.dt;
    otherwise
        error('displacement type unclear. Please check spelling')
end
  
  % Velocities from displacements
    VEL.U = DISP.Dx / time_increment;
    VEL.V = DISP.Dy / time_increment;
    VEL.W = DISP.Dz / time_increment;
    
 % total velocities in 2d and 3d
    VEL.V2d = sqrt(VEL.U.^2 + VEL.V.^2);
    VEL.V3d = sqrt(VEL.U.^2 + VEL.V.^2 + VEL.W.^2);

   [dUdy,dUdx] = gradient(VEL.U,COORDS.dy*GRID.ny,COORDS.dx*GRID.nx); % calculate dU/dx & dU/dy
   [dVdy,dVdx] = gradient(VEL.V,COORDS.dy*GRID.ny,COORDS.dy*GRID.ny); % calculate dV/dx & dV/dy

  % Strain-rate components
    VEL_DEVS.exx = dUdx;
    VEL_DEVS.exy = 0.5*(dUdy+dVdx);
    VEL_DEVS.eyx = 0.5*(dVdx+dUdy);
    VEL_DEVS.eyy = dVdy;

  % Vorticity
    VEL_DEVS.omega = 1/2*(dVdx - dUdy);
    
  % Vector angles
    [VEL_DEVS.theta, ~] = cart2pol(VEL.U,VEL.V);
    
  % Angular velocity
    VEL_DEVS.ang_vel = VEL_DEVS.theta / time_increment;
    
  % Invariants
    VEL_DEVS.i1 = VEL_DEVS.exx  +  VEL_DEVS.eyy;                                % tr(S)
    VEL_DEVS.i2 = VEL_DEVS.exx .*  VEL_DEVS.eyy - VEL_DEVS.exy .* VEL_DEVS.eyx; % det(S)
    
  % Angle of principal stretching axes
    VEL_DEVS.theta_p = 0.5*atan(2*VEL_DEVS.exy ./ (VEL_DEVS.exx - VEL_DEVS.eyy));
    
  % Principal stretching axes
    A = VEL_DEVS.exx + VEL_DEVS.eyy;
    B = VEL_DEVS.exx - VEL_DEVS.eyy;
    C = VEL_DEVS.exy;
    
    VEL_DEVS.emin = 0.5*(A) - ((0.5*(B)).^2 + C.^2).^0.5;
    VEL_DEVS.emax = 0.5*(A) + ((0.5*(B)).^2 + C.^2).^0.5;
    
  % max shear strain
    VEL_DEVS.gmax = VEL_DEVS.emax - VEL_DEVS.emin;
end

