% EURODYN FIGURES
% load('FEM_sol_h05_r_h2');
load('FEM_sol_h182_r_h2');
close all;
t=[0:ddt:(inData.Nper)*T]; %[sec]
heave = inData.a3.*sin(inData.omega3.*t+inData.phase3);

% figure;hold on;grid on;
% plot(t,heave, 'k')

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

figure;
% subplot(1,2,1);
hold on;grid on;
for ii = 1:length(t)
    w = solution.w(:,ii);
    sol(ii) = w(BBnodesTIP(NODE));
%     plot(t(ii)/inData.T3,heave,'r^');
end
% plot(t/inData.T3,sol/inData.a3,'b.-');
h1=plot(t/inData.T3,(sol+heave)/inData.a3,'k--');
% plot(t/inData.T3,heave/inData.a3,'r');
xlabel('$$time\,(sec)$$', 'interpreter','latex');
ylabel('$\alpha_{tip}/\alpha_{root}$', 'interpreter','latex');

% figure;
% % subplot(1,2,2);hold on;
% % plot3(pp(1,BBnodesTIP),pp(2,BBnodesTIP),w(BBnodesTIP),'ks','MarkerSize',3);
% % plot3(pp(1,BBnodesTIP(NODE)),pp(2,BBnodesTIP(NODE)),w(BBnodesTIP(NODE)),'ro','MarkerSize',10);
% % hh=pdeplot(pp,ee,tt,'XYData',w,"ZData",w,'colormap','jet');
% % colorbar;shading interp;view([25 25]);%axis equal;

error('er')
load('FEM_sol_h05_r_h2');
% subplot(1,2,1);
hold on;grid on;
for ii = 1:length(t)
    w = solution.w(:,ii);
    sol(ii) = w(BBnodesTIP(NODE));
%     plot(t(ii)/inData.T3,heave,'r^');
end
% plot(t/inData.T3,sol,'b.-');
h2=plot(t/inData.T3,(sol+heave)/inData.a3,'k.-');
% plot(t/inData.T3,heave,'r');
% xlabel('t [s]');ylabel('w tip')

xlim([1 3])

load('root.txt');
h4=plot(root(:,1)+1.25,root(:,2),'^-', 'LineWidth',1.5, 'MarkerSize',2);

load('inflexible_tip.txt');
h5=plot(inflexible_tip(:,1)+1.25,inflexible_tip(:,2),'o', 'LineWidth',1.5, 'MarkerSize',2);
h55=plot(root(:,1)+1.25,interp1(inflexible_tip(:,1),inflexible_tip(:,2),root(:,1),'pchip'),'o-',...
    'LineWidth',1.5, 'MarkerSize',2,'Color',h5.Color);

load('flexible_tip.txt');
h6=plot(flexible_tip(:,1)+1.25,flexible_tip(:,2),'s', 'LineWidth',1.5, 'MarkerSize',2);
h66=plot(root(:,1)+1.25,interp1(flexible_tip(:,1),flexible_tip(:,2),root(:,1),'pchip'),'o-',...
'LineWidth',1.5, 'MarkerSize',2,'Color',h6.Color);

legend([h1 h2 h4 h55 h66],'flexible (BEM-FEM)','rigid (BEM)',...
    'root (exp)', 'inflexible (exp)','flexible (exp)', 'interpreter','latex');

set(gca,'FontSize',15);
