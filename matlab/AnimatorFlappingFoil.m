function [] = AnimatorFlappingFoil(FOIL, PARAMETERS, KINESIS...
                         ,CFX, CFY, CM, Cp, EvolvedWake, Name)
% if animationSwitch = true, it produces a visualization of the 
% flapping foil movemennt. Giving information about the forcing
% CL and the pressure-type Kutta evaluation at each time step.

t = 0; % Initializing // corresponding to d=1

%% Retrieving Data
Nb = FOIL.Nb;
c = FOIL.c;
%NACAcode = '0012';%dataInput.Geom.NACAcode;
XR = FOIL.XR;

Xm = FOIL.Collocation(1, :);

Nper = PARAMETERS.Time.Nper; 
T = PARAMETERS.Time.T; 
N = PARAMETERS.Time.N;
dt = PARAMETERS.Time.dt;
U = PARAMETERS.U;
w = PARAMETERS.wh;

hdivc = PARAMETERS.HDC;
Str = PARAMETERS.STR;
pitchAmplitude = PARAMETERS.theta0_deg;
phaseDifference = PARAMETERS.psi_deg;


Xbt = KINESIS.Xbt;
Ybt = KINESIS.Ybt;
time = KINESIS.time;

yup = 3.5; %PLOT
mmmm=1 ; 
PlottingParameters = [yup, mmmm];

%% animation settings      
yup = PlottingParameters(1);

hfig=figure;
set(hfig,'Position',[1 41 1366 652],'PaperPositionMode','auto');
hax = subplot(2,1,1,'parent',hfig,'XGrid','on','YGrid','on'); %It returns the objects identity
hold on
% axis(hax,'equal');
set(hax,'XLim',[-Nper*U*T-XR*c XR*c+c],'YLim',[-5*max(max(KINESIS.Ybt)) 5*max(max(KINESIS.Ybt))]);
xlabel(hax,'x'); ylabel(hax,'y');

htit = title(hax,['CASE:   Str=',num2str(Str, '%5.3f'),',   \eta =',num2str(hdivc),',   \theta =',num2str(pitchAmplitude),' deg,   \psi =',num2str(phaseDifference),' deg,   \omega = ', num2str(w,'%5.2f'),' rad/s,   N_{p}=',num2str(Nb),'   N_{per}=',num2str(Nper),'   ;/T=',num2str(t/T, '%5.3f')]);

scale = max([max(max(CFX)),max(max(CFY)), max(max(CM))]);
hax2=subplot(2,2,3,'parent',hfig,'XGrid','on','YGrid','on');
hold on;
set(hax2,'XLim',[0 Nper],'YLim',[-scale*1.5 scale*1.5])
xlabel(hax2,'t/T');
ylabel(hax2,'C_F');
hold on
htit2 = title(hax2,'Non-dimensional coefficients');

hax3=subplot(2,4,7,'parent',hfig,'XGrid','on','YGrid','on');
hold on;
set(hax3,'XLim',[0 c],'YLim',[-max(max(Cp)) max(max(Cp))])
xlabel(hax3,'x');
ylabel(hax3,'C_p');
htit3 = title(hax3,['t/T=',num2str(t/T, '%5.3f')]);
hold on

hax4=subplot(2,4,8,'parent',hfig,'XGrid','on','YGrid','on');
hold on;
set(hax4,'XLim',[0 Nper],'YLim',[-0.1 0.1])
xlabel(hax4,'t/T');
ylabel(hax4,'');
htit4 = title(hax4,['Pressure-type Kutta: t/T=',num2str(t/T, '%5.3f')]);

hold on

writerObj = VideoWriter([Name,'_Str',num2str(Str,'%5.2f'),'_hdc',num2str(hdivc),'_theta0',num2str(pitchAmplitude),'_Nb',num2str(Nb),'_Nper',num2str(Nper),'.avi'],'Motion JPEG AVI');   
open(writerObj);

%% lines
hbody=line(Xbt(1,:),Ybt(1,:),'Parent',hax,'Color','g','LineStyle','-');
hbfill=fill(Xbt(1,:),Ybt(1,:),'g','Parent',hax);

hw = line(Xbt(1,1),Ybt(1,1),'Parent',hax,'Color','k','LineStyle','--'); % a point

hCL=line(time(1:(1)),CFY(1:(1)),'Parent',hax2,'Color','r','LineStyle','-');
%     hCT_8=line(time(1:(1)),-CFX_8(1:(1)),'Parent',hax2,'Color','r','LineStyle','--');
hCT=line(time(1:(1)),CFX(1:(1)),'Parent',hax2,'Color','r','LineStyle','--');
hCM=line(time(1:(1)),CM(1:(1)),'Parent',hax2,'Color','r','LineStyle',':');
legend(hax2,'C_L','C_T','C_M','Location','NorthWest');

hCplower=line((Xm(Nb/2:-1:1)+XR*c),Cp(1,Nb/2:-1:1),'Parent',hax3,'Color','r','LineStyle','-');
hCpupper=line((Xm(Nb/2+1:1:Nb)+XR*c),Cp(1,Nb/2+1:1:Nb),'Parent',hax3,'Color','b','LineStyle','-');

hCpfem=line((Xm(Nb/2+1:1:Nb)+XR*c),Cp(1,Nb/2:-1:1)-Cp(1,Nb/2+1:1:Nb),'Parent',hax3,'Color','m','LineStyle','-');

legend(hax3,'Cp LowerSide','Cp UpperSide','FEM forcing');

hdCp=line(time(1:(1)),(Cp(1:(1),Nb)-Cp(1:(1),1))','Parent',hax4,'Color','r','LineStyle','-');
    

for d = 2: N+1
    %% Since I want to reproduce the final results here, I need to make
    % sure that in each time step I reveal the right amount of information
    % based on what was known until then.
    % Check out(!): ExistingUntilThatTime

    t=(d-1)*dt;

    set(hbody,'XData',Xbt(d,:),'YData',Ybt(d,:));
    set(hbfill,'XData',Xbt(d,:),'YData',Ybt(d,:));
 
    set(hCL,'XData',time(1:(d)),'YData',CFY(1:(d)));
    set(hCT,'XData',time(1:(d)),'YData',CFX(1:(d)));
    set(hCM,'XData',time(1:(d)),'YData',CM(1:(d)));
    
    set(hCplower,'YData',Cp(d,Nb/2:-1:1),'Color','r');
    set(hCpupper,'YData',Cp(d,Nb/2+1:1:Nb),'Color','b');
    set(hCpfem,'YData',Cp(d,Nb/2:-1:1)-Cp(d,Nb/2+1:1:Nb),'Color','k');


    set(hdCp,'XData',time(1:d),'YData',(Cp((1:(d)),Nb)-Cp((1:(d)),1))');

    
    htit = title(hax,['CASE:    Str=',num2str(Str, '%5.3f'),',   \eta =',num2str(hdivc),',   \theta =',num2str(pitchAmplitude),' deg,   \psi =',num2str(phaseDifference),' deg,   \omega = ', num2str(w,'%5.2f'),' rad/s,   N_{p}=',num2str(Nb),'   N_{per}=',num2str(Nper),'   t/T=',num2str(t/T, '%5.3f')]);
    htit3 = title(hax3,['t/T=',num2str(t/T, '%5.3f')]);
    htit2 = title(hax2,'Non-dimensional coefficients');
    htit4 = title(hax4,['\it Pressure-type Kutta: t/T=',num2str(t/T, '%5.3f')]);

    %% --------------------------------------------------------------------
    ExistingUntilThatTime = (d-1) - 1; %(d-1) - 1; %Kutta strip does not count.
    EvolvedWake.Nf = ExistingUntilThatTime;
    %% --------------------------------------------------------------------
    
    Xfplot = flip(EvolvedWake.Xf(end:-1:end-EvolvedWake.Nf));
    Yfplot = flip(EvolvedWake.Yf(end:-1:end-EvolvedWake.Nf));
    mfplot = flip(EvolvedWake.mf(end:-1:end-EvolvedWake.Nf));
    % Actual size known 1: (Nf_8) new ones plus Kutta strip

    Xfp(1)=Xbt(d,1);
    Xfp(2:EvolvedWake.Nf+2)= Xfplot;
    
    % Actual size known 1: (Nf_8 + 1) new ones plus Kutta strip
    Yfp(1)=Ybt(d,1);
    Yfp(2:EvolvedWake.Nf+2)=Yfplot;
    
    % Actual size known 1: (Nf_8 + 1) new ones plus Kutta strip
    set(hw,'XData',Xfp(1:EvolvedWake.Nf+2),'YData',Yfp(1:EvolvedWake.Nf+2));

    Xfm = (Xfp(1:(EvolvedWake.Nf+1))+Xfp(2:(EvolvedWake.Nf+2)))/2;
    Yfm = (Yfp(1:(EvolvedWake.Nf+1))+Yfp(2:(EvolvedWake.Nf+2)))/2;
    
    DXf = Xfp(2:(EvolvedWake.Nf+2))-Xfp(1:(EvolvedWake.Nf+1));
    DYf = Yfp(2:(EvolvedWake.Nf+2))-Yfp(1:(EvolvedWake.Nf+1));
    
    DISTf=(DXf.^2+DYf.^2).^(1/2);
    SINTHEf=DYf./DISTf;
    COSTHEf=DXf./DISTf;

    nxf = -mfplot.*SINTHEf;
    nyf = mfplot.*COSTHEf;

    mmmm = PlottingParameters(2);
    if d == 2
        hnw = quiver(Xfm(1:EvolvedWake.Nf+1),Yfm(1:EvolvedWake.Nf+1),mmmm*nxf(1:EvolvedWake.Nf+1),mmmm*nyf(1:EvolvedWake.Nf+1),0,'Color','r','Parent',hax);     
    else
        set(hnw,'XData',Xfm(1:EvolvedWake.Nf+1),'YData',Yfm(1:EvolvedWake.Nf+1),'UData',mmmm*nxf(1:EvolvedWake.Nf+1),'VData',mmmm*nyf(1:EvolvedWake.Nf+1));
    end


    Fig=getframe(hfig);
    writeVideo(writerObj,Fig);
end

close(writerObj);
end

