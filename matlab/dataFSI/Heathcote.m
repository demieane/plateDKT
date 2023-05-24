%% POST-PROCESSING;

% u_fromC=Usol(1:GEN_fromC); % the vector of nodal unknowns (w1;bx1;by1;....wN;bxN;byN)
% 
% w=u_fromC(1:3:end);   % vertical displacement
% bx_fromC=u_fromC(2:3:end);  % rotation x
% by_fromC=u_fromC(3:3:end);  % rotation y

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

% load newmark;
% load implicitEuler;
% load crankNicolson;

% load FEM_sol_h182_r_h2_NEW;
% load solFromC;

figure;
hold on;grid on;
for ii = 2:length(t)
% %     wNewmark = solutionNewmark.w(:,ii);
% %     h1=plot(t(ii)/inData.T3,wNewmark(BBnodesTIP(NODE)),'bo', 'MarkerSize',2);
% %     %
% %     wImplicitEuler = solutionImplicitEuler.w(:,ii);
% %     h2=plot(t(ii)/inData.T3,wImplicitEuler(BBnodesTIP(NODE)),'r^', 'MarkerSize',2);
% %     %
% %     wCrankNicolson = solutionCrankNicolson.w(:,ii);
% %     h3=plot(t(ii)/inData.T3,wCrankNicolson(BBnodesTIP(NODE)),'ks', 'MarkerSize',2);
    %
    w = solution.w(:,ii);
    h0=plot(t(ii)/inData.T3,w(BBnodesTIP(NODE)),'go-', 'MarkerSize',2);
    %
%     w_fromC = w;%-w99(:,ii);
%     h99=plot(t(ii)/inData.T3,w_fromC(BBnodesTIP(NODE)),'m^-', 'MarkerSize',2);
end
legend([h0 ],'solution')
% legend([h0 h1 h2 h3 h99],'solution', 'newmark', 'implicit euler','crank-nicolson','c')
xlabel('t [s]');ylabel('w tip');


figure;
hold on; grid on;
plot3(pp(1,BBnodesTIP),pp(2,BBnodesTIP),w(BBnodesTIP),'ks','MarkerSize',3);
plot3(pp(1,BBnodesTIP(NODE)),pp(2,BBnodesTIP(NODE)),w(BBnodesTIP(NODE)),'ro','MarkerSize',10);
hh=pdeplot(pp,ee,tt,'XYData',w,"ZData",w,'colormap',viridis);
colorbar;shading interp;view([25 25]);%axis equal;


% C ~ Mglob, Kglob.
% In case load disappears the displacement needs to go to zero.
