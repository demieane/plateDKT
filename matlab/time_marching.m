%%% one-degree oscillator
% q'' + wo^2 q = 0
%
wo = pi; %rad/s




% set u = [q', q] to reduce the order of the system to 
% u' = A u 
A = [0 -wo^2;  1    0];

%analytic solution
T = 3;%s
h = T/32;
t = [0:h:T];
q = cos(wo*t);
% -wo^2 cos(wo t) + wo^2 * cos(wt) = 0
u=zeros(2,length(t));
u(:,1) = [0;1]; %[q'(0), q(0)]

figure;
h1=plot(t,q,'--')

for ii = 1:length(t)
    f(:,ii) = A*u(:,ii); %f(t,y(t))
    u(:,ii+1)= u(:,ii) + h*f(:,ii);%forward euler
    u_f_euler(:,ii) = u(:,ii+1);

    fguess = A*u(:,ii+1);
    u_b_euler(:,ii+1) = u(:,ii) + h*fguess;%backward euler;

    u(:,ii+1) = u(:,ii) + 0.5*h*(fguess + f(:,ii));
end
hold on; grid on;
h2=plot(t,u_f_euler(2,1:length(t)));
h3=plot(t,u_b_euler(2,1:length(t)));
h4=plot(t,u(2,1:length(t)), 'b.-');


M = 1;C=0;K=wo^2;
P = zeros(1,length(t));

gamma = 0.5;
beta = 0.25;
q=zeros(1,length(t));
q(1)=1;
qdot=zeros(1,length(t));
qdot(1)=0;
qdot2= zeros(1,length(t));
qdot2(1)=(u(1,2)-u(1,1))/h;
zeros(1,length(t));

for ii = 1:length(t)
    qdot(ii+1) = qdot(ii)+(1-gamma)*h*qdot2(ii);
    q(ii+1) = q(ii) + h*qdot(ii)+h^2*(1/2-beta)*qdot2(ii);
    qdot2(ii+1) = (P(ii) - C*qdot(ii+1) - K*q(ii+1))/(M-gamma*h*C - h^2*beta*K);
    % re-calculate
    qdot(ii+1) = qdot(ii)+(1-gamma)*h*qdot2(ii) + gamma*h*qdot2(ii+1);
    q(ii+1) = q(ii) + h*qdot(ii)+h^2*(1/2-beta)*qdot2(ii) + h^2*beta*qdot2(ii+1);
end


h5=plot(t,q(1,1:length(t)),'k--');

legend([h1 h2 h3 h4 h5], 'exact solution', 'forward euler',...
    'backward euler', 'crank-nicolson', 'newmark');
xlabel('time [sec]');
ylabel('response')