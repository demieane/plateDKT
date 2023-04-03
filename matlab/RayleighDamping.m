function [ C , res_Freq, a, b] = RayleighDamping( dens, xthvec, thicvec, Dvec, L, Kglob, Mglob, RayleighVisual)
%% ------------------------------------------------------------------------
% References : Theory of Vibrations McLachlan 1951 DOVER
% from NTUA library p.139
% -------------------------------------------------------------------------
% Calculate Eigenvalues
% Damping Matrix C, based on Rayleigh Method

% [XX,lamM,flag]=eigs(Kglob,Mglob,15,'sm');
% cc=sort(diag(lamM));
% res_Freq=sqrt(sort(diag(lamM),'ascend'))./(2*pi); % frequency :[1/sec]

[XX,lamM,flag]=eigs(Kglob,Mglob,15,'smallestabs',...
    'IsSymmetricDefinite',0, 'Tolerance',1^(-3),'MaxIterations',2);
cc=sort(diag(lamM))
res_Freq=sqrt(sort(diag(lamM),'ascend'))./(2*pi);


%% Eigenvalues - For modal reduction kratame autes pou antistoixoyn se mikroteres
% Dav = trapezoidal(xthvec,Dvec);
% thickav = trapezoidal(xthvec,thicvec);
% [bn] = eig(Kglob,Mglob); %Ascending order
% wn = sqrt(bn); %Natural/Resonance frequency [rad/s] in vacuuo
% Area = (2*thickav)*L;
% Fundamental_kl_squared = sqrt(bn*dens*Area*L^4/Dav); %(lk)^2 Fundamental frequency: Non-dimensional
% res_Freq = Fundamental_kl_squared/(2*pi); % frequency :[1/sec]

%--------------------------------------------------------------------------
% Damping: 0.1 at first mode, 0.22 at second mode
% Arguement #1: vector of natural frequencies, eigenfrequencies,resonance
% frequencies, whatever
% Arguement #2: m-th mode(here, m = 2)
% Arguement #3: n-th mode = 2.5m-th mode (here, n = 3)
% Arguement #4: Damping at first mode
% Arguement #5: Damping at m-th mode
% % % [a,b] = rayleighDampCoef(resFreq,2,3,1,0.5);
%--------------------------------------------------------------------------

% PARAZ_FREQUENCIES = 0;
% 
% if PARAZ_FREQUENCIES
%     mth_mode = 2;
%     nth_mode = 3;
%     % Ray4
% % %     zeta_1 = 0.020; 
% % %     zeta_m = 0.025; 
%     % Ray5
% % %     zeta_1 = 0.30; 
% % %     zeta_m = 0.27; 
% %*************BEST OPTIONS************************ 
%     % Ray3
%     zeta_1 = 0.3; 
%     zeta_m = 0.25; 
%     % Ray6
% %     zeta_1 = 0.29; 
% %     zeta_m = 0.27; 
% %*************************************************
%     RayleighFreqInput = [4.71, 6.22, 8.17];%wn;
%     [a,b] = rayleighDampCoef(RayleighFreqInput,mth_mode,nth_mode,zeta_1,zeta_m);
% else
    % BASED ON IN-VACUO EIGEN FREQUENCIES
    mth_mode = 2;
    NumberSignificantModes = 3;
    nth_mode = 5;%floor(2.5*NumberSignificantModes);
%     zeta_1 = 0.05;
    zeta_1 = 0.01;%0.1;
%     zeta_1divm = 2.5;
    zeta_m = 0.01;%0.1;%0.001;%zeta_1/zeta_1divm
    
%     zeta_1 = 0.020; 
%     zeta_m = 0.025; 
    [a,b] = rayleighDampCoef(res_Freq(1:end),mth_mode,nth_mode,zeta_1,zeta_m);
    % For constant thickenss Paraz et al. 2014 EXP a = 2.24, b=0.01
    res_Freq(1:nth_mode);
% end

if b<0
    b = b*(-1);
end

C = a*Mglob+b*Kglob;

% Plot zeta
n = nth_mode;

if RayleighVisual
    figure;
    omegaInput = res_Freq(1:n);
    omega_test = linspace(res_Freq(1),res_Freq(n),100);
    zetaInput = (a/2)./omegaInput + (b/2)*omegaInput;
    zeta_curve = (a/2)./omega_test + (b/2)*omega_test;
    
    plot(omega_test,zeta_curve,'-*'); grid on; hold on;
    
    plot(omegaInput,zetaInput,'*r'); hold on;
    
%     lala=find(omega_test>=RayleighFreqInput(n-1)); %To lala einai mia endiamesh timh pou elegxw an paremballetai swsta
%     plot(RayleighFreqInput(n-1),zeta_curve(lala(1)),'*r'); hold on;
%     plot(RayleighFreqInput(n),zeta_curve(end),'*r'); hold on;
    title('Rayleigh Damping')
end

% res_Freq(1:5)
a
b
% error('er')
% if RayleighVisual
%     figure;
%     omega_test = linspace(res_Freq(1),res_Freq(n),100);
%     
%     zeta_curve = (a/2)./omega_test + (b/2)*omega_test;
%     
%     plot(omega_test,zeta_curve,'-*'); grid on; hold on;
%     
%     plot(res_Freq(1),zeta_curve(1),'*r'); hold on;
%     
%     lala=find(omega_test>=res_Freq(n-1)); %To lala einai mia endiamesh timh pou elegxw an paremballetai swsta
%     plot(res_Freq(n-1),zeta_curve(lala(1)),'*r'); hold on;
%     plot(res_Freq(n),zeta_curve(end),'*r'); hold on;
%     title('Rayleigh Damping')
% end

% fprintf(' DAMPING PARAMETERS:\n',L)
% % fprintf(' -1rst Resonant Frequency: res_Freq= %5.3f Hz \n',res_Freq(1))
% fprintf(' -1rst Angular Natural Frequency: wmega_natural= %5.3f rad/sec \n',wn(1))
% fprintf(' -1rst Angular Fundamental Frequency: (kl)^2= %5.3f  \n',Fundamental_kl_squared(1))
% fprintf(' -Rayleigh Parameter1: a = %5.3f  \n',a)
% fprintf(' -Rayleigh Parameter2: b = %5.3f  \n\n',b)

end

