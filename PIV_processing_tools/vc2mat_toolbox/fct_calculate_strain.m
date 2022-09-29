function PLT = fct_calculate_strain(INPUT,nx,ny,F,whole_mask)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

switch INPUT.strain_type
    case 'infinitesimal'
        e = NaN(2,2,nx,ny);
        
        for ix = 1:nx
            for iy = 1:ny
                if ~isnan(whole_mask(ix,iy))
                    
                    % get local deformation gradient tensor
                    F_loc = F(:,:,ix,iy);
                    
                    % calculate local infinitesimal strain components
                    e_loc = 0.5 * (F_loc + F_loc') - eye(2,2);
                    e(:,:,ix,iy) = e_loc;
                    
                end
            end
        end
        
        PLT.vorticity = squeeze(0.5 * (F(2,1,:,:) - F(1,2,:,:)));
        PLT.exx       = squeeze(e(1,1,:,:));
        PLT.exy       = squeeze(e(1,2,:,:));
        PLT.eyx       = squeeze(e(2,1,:,:));
        PLT.eyy       = squeeze(e(2,2,:,:));
        
    case 'finite'
        E = NaN(2,2,nx,ny);
        
        for ix = 1:nx
            for iy = 1:ny
                if ~isnan(whole_mask(ix,iy))
                    
                    % get local deformation gradient tensor
                    F_loc = F(:,:,ix,iy);
                    
                    % calculate local finite strain components
                    E_loc = 0.5 * (F_loc' * F_loc - eye(2,2));
                    E(:,:,ix,iy) = E_loc;
                    
                end
            end
        end
        
        PLT.Exx       = squeeze(E(1,1,:,:));
        PLT.Exy       = squeeze(E(1,2,:,:));
        PLT.Eyx       = squeeze(E(2,1,:,:));
        PLT.Eyy       = squeeze(E(2,2,:,:));
        
    case 'stretch'
        U = NaN(2,2,nx,ny);
        V = NaN(2,2,nx,ny);
        R = NaN(2,2,nx,ny);
        A = NaN(nx,ny);
        
        for ix = 1:nx
            for iy = 1:ny
                if ~isnan(whole_mask(ix,iy))
                    
                    % get local deformation gradient tensor
                    F_loc = F(:,:,ix,iy);
                    
                    % calculate stretches
                    [R_loc, U_loc, V_loc, ~, RotAngleDeg_now] = fct_PolarDecomposition(F_loc);
                    
                    R(:,:,ix,iy) = R_loc;
                    U(:,:,ix,iy) = U_loc;
                    V(:,:,ix,iy) = V_loc;
                    A(ix,iy)     = RotAngleDeg_now;
                end
            end
        end
        
        PLT.U11 = squeeze(U(1,1,:,:));
        PLT.U12 = squeeze(U(1,2,:,:));
        PLT.U21 = squeeze(U(2,1,:,:));
        PLT.U22 = squeeze(U(2,2,:,:));
        
        PLT.V11 = squeeze(V(1,1,:,:));
        PLT.V12 = squeeze(V(1,2,:,:));
        PLT.V21 = squeeze(V(2,1,:,:));
        PLT.V22 = squeeze(V(2,2,:,:));
        
        PLT.Ang = A;

    otherwise
end

end

