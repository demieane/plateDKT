%% POST-PROCESSING;
load('solMatlab_dynamic_test_newmark');
wNewmark = solution.w;

load('solMatlab_dynamic_test_crankNicolson');
wCrankNicolson = solution.w;

load('solMatlab_dynamic_test_implicitEuler');
wImplicitEuler = solution.w;

NODE = 200;

figure;
hold on;grid on;
for ii = 2:length(t)
    wNM_data(ii) = wNewmark(NODE,ii);
    wCR_data(ii) = wCrankNicolson(NODE,ii);
    wIE_data(ii) = wImplicitEuler(NODE,ii);
%     h0=plot(t(ii)/inData.T3,w(ii),'b.-', 'MarkerSize',2);
end
h1 = plot(t/inData.T3,wNM_data,'b.-', 'MarkerSize',4);
hold on;grid on;
h2 = plot(t/inData.T3,wCR_data,'ro', 'MarkerSize',4);
h3 = plot(t/inData.T3,wIE_data,'ks', 'MarkerSize',4);
legend([h1 h2 h3],'newmark (MATLAB)', 'Crank-Nicolson (MATLAB)','Implicit Euler (MATLAB)');
xlabel('t/T');
ylabel(['w, NODE=',num2str(NODE)]);


error('er')


% u_fromC=Usol(1:GEN_fromC); % the vector of nodal unknowns (w1;bx1;by1;....wN;bxN;byN)
% 
% w=u_fromC(1:3:end);   % vertical displacement
% bx_fromC=u_fromC(2:3:end);  % rotation x
% by_fromC=u_fromC(3:3:end);  % rotation y

% load solution_newmark;
% % load implicitEuler;
% load solution_crankNicolson;
% load solution_TEST;

% load solution_BENCH_matlab

e = ee;
p = pp;
% solutionNewmark  = solution;
% save newmark solutionNewmark
% solutionImplicitEuler  = solution;
% save implicitEuler solutionImplicitEuler
% solutionCrankNicolson  = solution;
% save crankNicolson solutionCrankNicolson

% load FEM_sol2
% BnodesTIP =find(e(5,:)==3) 
BnodesTIP =find(e(5,:)==1) 
%
BBnodesTIP = BnodesTIP.*0;
for i=1:length(BnodesTIP)
    BBnodesTIP(i)=e(1,BnodesTIP(i));
end
BBnodesTIP=sort(BBnodesTIP);  
% Bdofs1TIP=ID(1,BBnodesTIP);%w
% Bdofs2=ID(2,BBnodesTIP);%bx
% Bdofs3=ID(3,BBnodesTIP);%by
NODE = 1;%round(length(BBnodesTIP));




% load FEM_sol_h182_r_h2_NEW;
% load solFromC;

figure;
hold on;grid on;
for ii = 2:length(t)
%     wNewmark = solutionNewmark.w(:,ii);
%     h1=plot(t(ii)/inData.T3,wNewmark(BBnodesTIP(NODE)),'bo', 'MarkerSize',2);
% %     %
% %     wImplicitEuler = solutionImplicitEuler.w(:,ii);
% %     h2=plot(t(ii)/inData.T3,wImplicitEuler(BBnodesTIP(NODE)),'r^', 'MarkerSize',2);
% %     %
%     wCrankNicolson = solutionCrankNicolson.w(:,ii);
%     h3=plot(t(ii)/inData.T3,wCrankNicolson(BBnodesTIP(NODE)),'ks', 'MarkerSize',2);
    %
    w = solution.w(:,ii);
    h0=plot(t(ii)/inData.T3,w(BBnodesTIP(NODE)),'go-', 'MarkerSize',2);
    %
    w_fromC = 0;%solutionC.w(:,ii);
    h99=plot(t(ii)/inData.T3,w_fromC(BBnodesTIP(NODE)),'m^-', 'MarkerSize',2);
end
legend([h0 h99],'newmark (MATLAB)','newmark (C)');
% legend([h0 h1 h2 h3 h99],'solution', 'newmark', 'implicit euler','crank-nicolson','c')
xlabel('t [s]');ylabel('w tip');


figure;
subplot(1,2,1);
hold on; grid on;
plot3(pp(1,BBnodesTIP),pp(2,BBnodesTIP),w(BBnodesTIP),'ks','MarkerSize',3);
plot3(pp(1,BBnodesTIP(NODE)),pp(2,BBnodesTIP(NODE)),w(BBnodesTIP(NODE)),'ro','MarkerSize',10);
hh=pdeplot(pp,ee,tt,'XYData',w,"ZData",w,'colormap',viridis);
colorbar;shading interp;view([25 25]);
%axis equal;
title('matlab');

% figure;
subplot(1,2,2);
hold on; grid on;
plot3(pp(1,BBnodesTIP),pp(2,BBnodesTIP),w_fromC(BBnodesTIP),'ks','MarkerSize',3);
plot3(pp(1,BBnodesTIP(NODE)),pp(2,BBnodesTIP(NODE)),w_fromC(BBnodesTIP(NODE)),'ro','MarkerSize',10);
hh=pdeplot(pp,ee,tt,'XYData',w_fromC,"ZData",w_fromC,'colormap',viridis);
colorbar;shading interp;view([25 25]);%axis equal;
title('c');


% C ~ Mglob, Kglob.
% In case load disappears the displacement needs to go to zero.
