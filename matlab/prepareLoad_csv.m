function [] = prepareLoad_csv(x,y,fx,filename,ID)
%==========================================================================
% Preparation of .csv file for External Data ANSYS
%==========================================================================
N1=size(x,1);
N2=size(x,2);

z=0;
A=zeros(N1*N2,3);
for ii=1:N1
    for jj=1:N2
        z=z+1;
        A(z,1)=x(ii,jj);
        A(z,2)=y(ii,jj);
        A(z,3)=fx(ii,jj);
    end
end
% writecell({'x(m)','y(m)','pressure(Pa)'},'load.csv');
% writematrix(A,'load.csv','WriteMode','append');
if ID ==1
    writecell({'x(m)','y(m)','pressure(Pa)'},filename);
elseif ID ==2
     writecell({'x(m)','y(m)','thickness(m)'},filename);
end
writematrix(A,filename,'WriteMode','append');
end

