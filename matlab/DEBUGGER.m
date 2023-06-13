
% % load mat
% % 
% % max(max(abs(A - Adebug)))
% % max(max(abs(B - Bdebug)))
% % 
% % max(max(abs(B(1:sizeM,1:sizeM) - Bdebug(1:sizeM,1:sizeM))))
% % max(max(abs(B(1:sizeM,sizeM+1:end) - Bdebug(1:sizeM,sizeM+1:end))))
% % 
% % max(max(abs(full(Kglob) - -Bdebug(1:sizeM,sizeM+1:end))))

% file = fopen('test_lin_solve', 'wb');
% % AA
% fwrite(file, size(AA,1),'int');
% fwrite(file, size(AA,2),'int');
% for ii = 1:size(AA,1)
%     for jj = 1:size(AA,2)
%         fwrite(file, AA(ii,jj),precision);
%     end
% end
% % rhs
% fwrite(file, size(rhs,1),'int');
% fwrite(file, size(rhs,2),'int');
% for ii = 1:size(rhs,1)
%     for jj = 1:size(rhs,2)
%         fwrite(file, rhs(ii,jj),precision);
%     end
% end
% fclose(file);



%%
% Read solution from binary file
modeFem = 1;
precision = 'double';

% fileID = fopen('../c/OUTDATA_FEM_DEBUG.bin','rb')
fileID = fopen('../c/OUTDATA_FEM_double.bin','rb')
GEN_fromC = fread(fileID,1,'int')
rowsUsol = fread(fileID,1,'int')
colsUsol = fread(fileID,1,'int')

for i = 1:rowsUsol
    for j = 1:colsUsol
%         Cdamp(i,j)=fread(fileID,1,precision);
%         Gdebug(i,j)=fread(fileID,1,precision);
%         AAdebug(i,j)=fread(fileID,1,precision);
%         qdot2_debug(i,j)=fread(fileID,1,precision);
        Usol(i,j)=fread(fileID,1,precision);
    end
end

fclose(fileID);

solutionC = struct( 'w',[], 'bx',[], 'by', [],...
        'w_dot',[], 'bx_dot',[], 'by_dot', [],...
        'uu', [], 'uu_dot',[]); 

load solution_BENCH_matlab

    
for d = 1:length(t)-1
    [solutionC] = solutionRetriever(GEN_fromC, sizeM, d+1, length(t), Usol, solutionC);%[w,bx,by]
end


% % rowsUsol = fread(fileID,1,'int')
% % colsUsol = fread(fileID,1,'int')
% % 
% % for i = 1:rowsUsol
% %     for j = 1:colsUsol
% % %         Cdamp(i,j)=fread(fileID,1,precision);
% % %         Gdebug(i,j)=fread(fileID,1,precision);
% %         rhs(i,j)=fread(fileID,1,precision);
% %     end
% % end
% % 
% % 
% % rowsUsol = fread(fileID,1,'int')
% % colsUsol = fread(fileID,1,'int')
% % 
% % for i = 1:rowsUsol
% %     for j = 1:colsUsol
% % %         Cdamp(i,j)=fread(fileID,1,precision);
% % %         Gdebug(i,j)=fread(fileID,1,precision);
% %         qdot2_buffer(i,j)=fread(fileID,1,precision);
% %     end
% % end

% % fclose(fileID);

% solutionC = struct( 'w',[], 'bx',[], 'by', [],...
%         'w_dot',[], 'bx_dot',[], 'by_dot', [],...
%         'uu', [], 'uu_dot',[]); 
%         
% for d = 1:length(t)-1
%     [solutionC] = solutionRetriever(GEN, sizeM, d+1, length(t), Usol, solutionC);%[w,bx,by]
% end
        
        

Usol = AAdebug\rhs;

Usol(1:20)'
qdot2_debug(1:20,1)'
qdot2_buffer(1:20,1)'

% Usol(sizeM-3:sizeM+10)'
% u(sizeM-3:sizeM+10,2)'

load test 
max(max(abs(qdot2_buffer - qdot2(:,1))))



% max(max(abs(full(C) - Cdamp)))
% C(end,end)

% max(max(abs(full(G(:,2)) - Gdebug(:,2))))
% max(max(abs(full(G) - Gdebug)))

error('er')

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
