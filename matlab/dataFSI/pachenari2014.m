function [txx] = pachenari2014(h, alpha, chord, span, x, y, IEN, p, e, t, debugOn)

N1=20;
N2=20;
[xx,yy]=meshgrid(linspace(0,chord,N1),linspace(0,span,N2));

h1 = h;
% r = 1;
% l = 0;
% alpha = 0.75;

x1 = 0;
x2 = chord;
y1 = h1*alpha;
y2 = h1;

slope = (y1-y2)/(x1-x2);
ct = y1-slope*x1;

tx_linear = slope .*xx  + ct;

%%==================================

figure;
surf(xx,yy,tx_linear);

% error('oo')
% Triangle center
xm=(1/3)*(x(IEN(1,:))+x(IEN(2,:))+x(IEN(3,:)));    % x barycentric coordinate
ym=(1/3)*(y(IEN(1,:))+y(IEN(2,:))+y(IEN(3,:)));    % y barycentric coordinate

size(xm);
size(ym);

%txx=interp2(xx,yy,tx_linear,xm,ym,'spline');%is this OK?
txx=griddata(xx,yy,tx_linear,xm,ym,'linear');%is this OK?
        
if debugOn
   
    figure;hold on;grid on;
    plot3(xm,ym,txx,'o');
    title('NEW matlab interpolated expression for txx');
    xlabel('x-axis');
    ylabel('y-axis');view([-25 25])

    figure
    hold on;plot3(xm,ym,txx,'ks')
   % pdeplot(p,e,t,'xydata',txx,'zdata',txx,'colormap','copper')
    xlabel('x-axis');
    xlabel('y-axis');
    xlabel('z-axis');grid on;hold on;view([-25 25])
    title('thickness thick(x,y)')
end

end

