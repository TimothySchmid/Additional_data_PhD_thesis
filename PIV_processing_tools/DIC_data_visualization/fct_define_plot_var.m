function PLT = fct_define_plot_var(INPUT,DISP,DISP_DEVS,VEL,VEL_DEVS)

% fct_define_plot_var
%
% DEFINE PLOTTING VARIABLE
%
% Input:
%   - INPUT      --> Structure containing initial input parameters
%   - DIPS       --> Displacement values
%   - DISP_DEV   --> Displacement derived values
%   - VEL      --> Velocity components
%   - VEL_DEVS --> Velocity derivative values
%
% ======================================================================= %

switch INPUT.plot_val
    
    % displacements
    case 'Dx'
        PLT.val = DISP.Dx;
        PLT.lab = 'Displacement Dx';
    case 'Dy'
        PLT.val = DISP.Dy;
        PLT.lab = 'Displacement Dy';
    case 'Dz'
        PLT.val = DISP.Dz;
        PLT.lab = 'Displacement Dz';
    case 'D2d'
        PLT.val = DISP.D2d;
        PLT.lab = 'Total 2D displacement';
    case 'D3d'
        PLT.val = DISP.D3d;
        PLT.lab = 'Total 3D displacement';
        
        % displacement derived
    case 'Exx'
        PLT.val = DISP_DEVS.Exx;
        PLT.lab = 'Normal strain E_{xx}';
    case 'Exy'
        PLT.val = DISP_DEVS.Exy;
        PLT.lab = 'Shear strain E_{xy}';
    case 'Eyx'
        PLT.val = DISP_DEVS.Eyx;
        PLT.lab = 'Shear strain E_{yx}';
    case 'Eyy'
        PLT.val = DISP_DEVS.Eyy;
        PLT.lab = 'Normal strain E_{yy}';
    case 'Gavg'
        PLT.val = DISP_DEVS.Gavg;
        PLT.lab = 'Average shear strain';
    case 'Div'
        PLT.val = DISP_DEVS.Div;
        PLT.lab = '2D Divergence';
    case 'I1'
        PLT.val = DISP_DEVS.I1;
        PLT.lab = 'First invariant of strain';
    case 'I2'
        PLT.val = DISP_DEVS.I2;
        PLT.lab = 'Second invariant of strain';
    case 'Theta_P'
        PLT.val = DISP_DEVS.Theta_P;
        PLT.lab = 'Rotation angle';
    case 'Emin'
        PLT.val = DISP_DEVS.Emin;
        PLT.lab = 'Principal shortening axis magnitude';
    case 'Emax'
        PLT.val = DISP_DEVS.Emax;
        PLT.lab = 'Principal stretching axis magnitude';
    case 'Gmax'
        PLT.val = DISP_DEVS.Gmax;
        PLT.lab = 'Maximum shear strain magnitude';
        
        % velocities
    case 'U'
        PLT.val = VEL.U;
        PLT.lab = 'Velocity U';
    case 'V'
        PLT.val = VEL.V;
        PLT.lab = 'Velocity V';
    case 'W'
        PLT.val = VEL.W;
        PLT.lab = 'Velocity W';
    case 'V2d'
        PLT.val = VEL.V2d;
        PLT.lab = 'Total 2D velocity';
    case 'V3d'
        PLT.val = VEL.V3d;
        PLT.lab = 'Total 3D displacement';
        
        % velocity derived
    case 'exx'
        PLT.val = VEL_DEVS.exx;
        PLT.lab = 'Normal strain rate e_{xx}';
    case 'exy'
        PLT.val = VEL_DEVS.exy;
        PLT.lab = 'Shear strain rate e_{xy}';
    case 'eyx'
        PLT.val = VEL_DEVS.eyx;
        PLT.lab = 'Shear strain rate e_{yx}';
    case 'eyy'
        PLT.val = VEL_DEVS.eyy;
        PLT.lab = 'Normal strain rate E_{yy}';
    case 'omega'
        PLT.val = VEL_DEVS.omega;
        PLT.lab = 'Vorticity';
    case 'theta'
        PLT.val = VEL_DEVS.theta;
        PLT.lab = 'Rotation angle';
    case 'ang_vel'
        PLT.val = VEL_DEVS.ang_vel;
        PLT.lab = 'Angular velocity';
    case 'i1'
        PLT.val = VEL_DEVS.i1;
        PLT.lab = 'First invariant of strain rate';
    case 'i2'
        PLT.val = VEL_DEVS.i2;
        PLT.lab = 'Second invariant of strain rate';
    case 'theta_p'
        PLT.val = VEL_DEVS.theta_p;
        PLT.lab = 'Rotation angle';
    case 'emin'
        PLT.val = VEL_DEVS.emin;
        PLT.lab = 'Maximum strain rate magnitude';
    case 'emax'
        PLT.val = VEL_DEVS.emax;
        PLT.lab = 'Minimum strain rate magnitude';
    case 'gmax'
        PLT.val = VEL_DEVS.gmax;
        PLT.lab = 'Maximum shear strain rate magnitude';
        
    otherwise
        error('Plotting variable unknown. Please check spelling')
end
end

