% Read solution from binary file
modeFem = 1;
precision = 'double';

fileID = fopen('../c/OUTDATA_FEM_DEBUG.bin','rb')
rowsUsol = fread(fileID,1,'int')
colsUsol = fread(fileID,1,'int')

for i = 1:rowsUsol
    for j = 1:colsUsol
        AA(i,j)=fread(fileID,1,precision);
    end
end

rowsUsol1 = fread(fileID,1,'int')
colsUsol1 = fread(fileID,1,'int')

for i = 1:rowsUsol1
    for j = 1:colsUsol1
        utemp2(i,j)=fread(fileID,1,precision);
    end
end

rowsUsol2 = fread(fileID,1,'int')
colsUsol2 = fread(fileID,1,'int')

for i = 1:rowsUsol2
    for j = 1:colsUsol2
        usol(i,j)=fread(fileID,1,precision);
    end
end
% % 
% % 
usolDEBUG = mldivide(AA,utemp2);


usolDEBUG(1:10)'
usol(1:10)'
usol(sizeM+1)
usolDEBUG(sizeM+1)



max(abs(usol - usolDEBUG))

% % 
% % load('FEM_sol_h182_r_h2');

% % u(sizeM+1,3)
% % 
% % u(sizeM-10:sizeM,3)'
% % usolDEBUG(sizeM-10:sizeM)'
% % usol(sizeM-4:sizeM+10)'
