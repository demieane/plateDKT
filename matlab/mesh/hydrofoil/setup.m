function [Nodtot,X,Y] = setup(NEC,TAU,EPSMAX,NACA,PTMAX)

%%%%%%%% SET COORDINATES OF NODES ON BODY SURFACE %%%%%%%

Nlower=floor(NEC/2);
Nupper=Nlower;
NPOINTS=Nlower;
SIGN=-1;
NSTART=0;

for NSURF=1:2
   for N=1:NPOINTS
       
%        NPOINTS;
%        error('er')
       ficos=pi*(N-1)/NPOINTS;
%         [NSURF, N, ficos]
       if ((SIGN==-1)&&(ficos>pi/2))||((SIGN==1)&&(ficos<pi/2))  %((SIGN==-1)&&(N>NPOINTS/2))||((SIGN==1)&&(N<NPOINTS/2))
%            disp('cosine')
           Z=0.5*(1-cos(ficos)); % cosine distribution
       else
%            disp('isospace')
           Z=(N-1)/NPOINTS; % Isospace distribution
       end
       I=NSTART+N;
       [X(I),Y(I)]=body(Z,SIGN,TAU,EPSMAX,NACA,PTMAX);
%        [X(I), Y(I)]
   end
   
   NPOINTS=Nupper;
   SIGN=1;
   NSTART=Nlower;
end

Nodtot=Nlower+Nupper;
X(Nodtot+1)=X(1);
Y(Nodtot+1)=Y(1);

end

