function [R, U, V, RotAngleRad, RotAngleDeg] = fct_PolarDecomposition(F)
          % Compute symmetric right Cauchy-Green deformation tensor C
            C = F' * F;
            
          % Find Eigenvalues of C
            mu1 = 0.5 * (C(1,1) + C(2,2) + sqrt(4*C(1,2)*C(2,1) + (C(1,1) - C(2,2))^2) );
            mu2 = 0.5 * (C(1,1) + C(2,2) - sqrt(4*C(1,2)*C(2,1) + (C(1,1) - C(2,2))^2) );
            
          % Compute invariants of Cauchy-Green tensor C
            IC  = mu1 + mu2;
            IIC = mu1 * mu2;
            
          % Compute invariants of right stretch tensor U
            IU  = sqrt(IC + 2 * sqrt(IIC));
            IIU = sqrt(IIC);
            
          % Compute right stretch tensor U
            U = (C + IIU * eye(2)) / IU;
            
          % Compute inverse of right stretch tensor U
            Uinv = -IU / (IIU * (IIU + IC) + IIC) * (C - (IIU + IC) * eye(2));
            
          % Compute rotation tensor R
            R = F * Uinv;
            
          % Compute left stretch tensor V
            V = F * R';
            
          % Compute rotation angle using the four-quadrant inverse tangent
            RotAngleRad = atan2(R(2,1),R(1,1));
            RotAngleDeg = RotAngleRad * 180/pi;
end
