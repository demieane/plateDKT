%% POST-PROCESSING;

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

figure;
hold on;grid on;
for ii = 2:length(t)
    w = solution.w(:,ii);
    plot(t(ii)/inData.T3,w(BBnodesTIP(NODE)),'bo');
end
xlabel('t [s]');ylabel('w tip')

figure;hold on;
plot3(pp(1,BBnodesTIP),pp(2,BBnodesTIP),w(BBnodesTIP),'ks','MarkerSize',3);
plot3(pp(1,BBnodesTIP(NODE)),pp(2,BBnodesTIP(NODE)),w(BBnodesTIP(NODE)),'ro','MarkerSize',10);
hh=pdeplot(pp,ee,tt,'XYData',w,"ZData",w,'colormap','jet');
colorbar;shading interp;view([25 25]);%axis equal;


% C ~ Mglob, Kglob.
% In case load disappears the displacement needs to go to zero.
