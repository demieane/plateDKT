function [ THICK,CAMBER,BETA ] = NACA45( Z,TAU,EPSMAX,NACA,PTMAX )
%EVALUATE THICKNESS AND CAMBER FOR NACA 4- OR -5 DIGIT AIRFOIL

THICK=0;
if (Z>=10^(-10))
    THICK=5*TAU*(0.2969*sqrt(Z)-Z*(0.126+Z*(0.3537-Z*(0.2843-Z*0.1015))));
end
if(EPSMAX~=0)
    if (NACA<=9999)
        if (Z<=PTMAX)
            CAMBER=(EPSMAX/(PTMAX^2))*(2*PTMAX-Z)*Z;
            DCAMDX=2*(EPSMAX/(PTMAX^2))*(PTMAX-Z);
        else
            CAMBER=(EPSMAX/(1-PTMAX)^2)*(1+Z-2*PTMAX)*(1-Z);
            DCAMDX=2*(EPSMAX/(1-PTMAX)^2)*(PTMAX-Z);
        end
    else
        if(Z<=PTMAX)
            CAMBER=EPSMAX*(1-Z);
            DCAMDX=-EPSMAX;
        else
            W=Z/PTMAX;
            CAMBER=EPSMAX*W*((W-3)*W+3-PTMAX);
            DCAMDX=EPSMAX*3*W*(1-W)/PTMAX;
        end
    end
    BETA=atan(DCAMDX);
else
    CAMBER=0;
    BETA=0;
end

%CAMBER=CAMBER
%THICK=THICK
end

