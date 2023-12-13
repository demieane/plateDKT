%% Read solution from binary file

fileID = fopen('../c/OUTDATA_FEM_DEBUG.bin','rb')
precision = 'double';

% fileID = fopen('../c/OUTDATA_FEM.bin','rb')
rowsUsol = fread(fileID,1,'int')
colsUsol = fread(fileID,1,'int')

for i = 1:rowsUsol
    for j = 1:colsUsol
        AA_C(i,j)=fread(fileID,1,precision);
    end
end

rowsUsol = fread(fileID,1,'int')
colsUsol = fread(fileID,1,'int')

for i = 1:rowsUsol
%     for j = 1:colsUsol
        BB_C(i)=fread(fileID,1,precision);
%     end
end

fclose(fileID);

% sol1=AA_C\BB_C';
% sol1(1:10)