function uNEW = timeIntegration(u, d, GEN, Mglob, Kglob, C, G, dt, theta)
% 
%  A = [sparse(GEN,GEN) Mglob; speye(GEN,GEN) sparse(GEN,GEN)];
%  B = -[Kglob C; sparse(GEN,GEN) -speye(GEN,GEN)];

sizeM=size(Mglob,1);

A = [Mglob, sparse(sizeM,sizeM); sparse(sizeM,sizeM), speye(sizeM,sizeM)];
B = -[C, Kglob; -speye(sizeM,sizeM), sparse(sizeM,sizeM)];

%lamda (1: implicit Euler, 1/2: Crank - Nicolson)

% AA =  A + theta*dt*B;
% BB =  A - (1 - theta)*dt*B;

AA =  A - theta*dt*B;
BB =  A + (1 - theta)*dt*B;

Q = (1 - theta)*dt*G(:,d-1) + (theta)*dt*G(:,d);

% Q(1:10)'
% Q(sizeM+1:sizeM+10)'
% 
% rcond(full(AA));

% disp('Using sparse matrices');
% t11=mldivide(AA,(BB*u(:,d-1) + Q));
% 
% t11(sizeM+1:sizeM+10)'
% disp('Using dense matrices');
% t22=mldivide(full(AA),full((BB*u(:,d-1) + Q)));

% t22(sizeM+1:sizeM+10)'


% error('er')
u(:,d) = AA\(BB*u(:,d-1) + Q);

% rhs = (BB*u(:,d-1) + Q);

% % %     file = fopen('test_lin_solve.bin', 'wb');
% % %     % AA
% % %     size(full(AA),1)
% % %     size(full(AA),2)
% % %     fwrite(file, size(full(AA),1),'int');
% % %     fwrite(file, size(full(AA),2),'int');
% % %     for ii = 1:size(full(AA),1)
% % %         for jj = 1:size(full(AA),2)
% % %             fwrite(file, full(AA(ii,jj)),'double');
% % %         end
% % %     end
% % %     % rhs
% % %     size(rhs,1)
% % %     size(rhs,2)
% % %     fwrite(file, size(rhs,1),'int');
% % %     fwrite(file, size(rhs,2),'int');
% % %     for ii = 1:size(rhs,1)
% % %         for jj = 1:size(rhs,2)
% % %             fwrite(file, rhs(ii,jj),'double');
% % %         end
% % %     end
% % %     fclose(file);

% rhsDEBUG = 0;
% ii = 1;
% for jj = 1:sizeM*2
%     rhsDEBUG(ii) = rhsDEBUG(ii) + BB(jj,ii);%*u(jj,d-1);
% end
% rhsDEBUG
% 
% utemp2 = BB*u(:,d-1) + Q;
% utemp2(1:10)' 



uNEW = u(:,d);
% u(1:10,d)'
% u(end-10:1:end,d)'

% error('er')

%% debugging zone
% fileID = fopen('../c/OUTDATA_FEM_DEBUG.bin','rb')
% rowsUsol = fread(fileID,1,'int')
% colsUsol = fread(fileID,1,'int')
% 
% for i = 1:rowsUsol
%     for j = 1:colsUsol
%         Adebug(i,j)=fread(fileID,1,'double');
%     end
% end
% 
% rowsUsol1 = fread(fileID,1,'int')
% colsUsol1 = fread(fileID,1,'int')
% 
% for i = 1:rowsUsol1
%     for j = 1:colsUsol1
%         Bdebug(i,j)=fread(fileID,1,'double');
%     end
% end
% 
% rowsUsol2 = fread(fileID,1,'int')
% colsUsol2 = fread(fileID,1,'int')
% 
% for i = 1:rowsUsol2
%     for j = 1:colsUsol2
%         AAdebug(i,j)=fread(fileID,1,'double');
%     end
% end
% 
% rowsUsol3 = fread(fileID,1,'int')
% colsUsol3 = fread(fileID,1,'int')
% 
% for i = 1:rowsUsol3
%     for j = 1:colsUsol3
%         BBdebug(i,j)=fread(fileID,1,'double');
%     end
% end
% 

% save mat sizeM A B Adebug Bdebug
% max(max(abs(A - Adebug)))
% max(max(abs(B - Bdebug)))

% max(max(abs(AA - AAdebug)))
% max(max(abs(BB - BBdebug)))
% 
% fclose(fileID)
% 
% error('er')

%%




% 
% % a = max(max(abs(AA(1:sizeM,1:sizeM) - AAdebug(1:sizeM,1:sizeM)))) %a
% % b = max(max(abs(AA(1:sizeM,sizeM+1:end) - AAdebug(1:sizeM,sizeM+1:end))))%b
% % c = max(max(abs(AA(sizeM+1:end,1:sizeM) - AAdebug(sizeM+1:end,1:sizeM)))) %c
% % d = max(max(abs(AA(sizeM+1:end,sizeM+1:end) - AAdebug(sizeM+1:end,sizeM+1:end)))) %d
% % %
% % a = max(max(abs(BB(1:sizeM,1:sizeM) - BBdebug(1:sizeM,1:sizeM)))) %a
% % b = max(max(abs(BB(1:sizeM,sizeM+1:end) - BBdebug(1:sizeM,sizeM+1:end))))%b
% % c = max(max(abs(BB(sizeM+1:end,1:sizeM) - BBdebug(sizeM+1:end,1:sizeM)))) %c
% % d = max(max(abs(BB(sizeM+1:end,sizeM+1:end) - BBdebug(sizeM+1:end,sizeM+1:end)))) %d

% size(BB)
% size(BBdebug)
% size(u(:,d-1))
% 
% utemp1 = BB*u(:,d-1);
% utemp1debug = BBdebug*u(:,d-1);
% 
% utemp11 = utemp1(1:10)
% utemp12 = utemp1debug(1:10)
% 
% utemp2 = utemp1 + Q;
% utemp2_print = utemp2(1:10)

% uNEW = u(:,d);

% error('er')
% % % % 
% % % % 
% % % % % save mat A B Adebug Bdebug sizeM
% % % % % max(max(abs(A - Adebug)))
% % % % % max(max(abs(B - Bdebug)))
% % % % % 
% % % % % a = max(max(abs(A(1:sizeM,1:sizeM) - Adebug(1:sizeM,1:sizeM)))) %a
% % % % % b = max(max(abs(A(1:sizeM,sizeM+1:end) - Adebug(1:sizeM,sizeM+1:end))))%b
% % % % % c = max(max(abs(A(sizeM+1:end,1:sizeM) - Adebug(sizeM+1:end,1:sizeM)))) %c
% % % % % d = max(max(abs(A(sizeM+1:end,sizeM+1:end) - Adebug(sizeM+1:end,sizeM+1:end)))) %d
% % % % % %%
% % % % % 
% % % % % a = max(max(abs(B(1:sizeM,1:sizeM) - Bdebug(1:sizeM,1:sizeM)))) %a
% % % % % b = max(max(abs(B(1:sizeM,sizeM+1:end) - Bdebug(1:sizeM,sizeM+1:end))))%b
% % % % % c = max(max(abs(B(sizeM+1:end,1:sizeM) - Bdebug(sizeM+1:end,1:sizeM)))) %c
% % % % % d = max(max(abs(B(sizeM+1:end,sizeM+1:end) - Bdebug(sizeM+1:end,sizeM+1:end)))) %d
% % % % 
% % % % Q = (1 - theta)*dt*G(:,d-1) + (theta)*dt*G(:,d);
% % % % % u(1:10,d-1)';
% % % % % utemp1 = BB*u(:,d-1);
% % % % % utemp2 = BB*u(:,d-1) + Q;
% % % % 
% % % % 
% % % % % max(max(abs(Q - Qdebug)))
% % % % % max(max(abs(utemp2 - utemp2debug)))
% % % % % 
% % % % % a = max(abs(Q(1:sizeM) - Qdebug(1:sizeM))) %a
% % % % % b = max(abs(utemp2(sizeM+1:end) - utemp2debug(sizeM+1:end))) %a%b
% % % % 
% % % % % error('er')
% % % % 
% % % % % utemp2=utemp2(783:783+10)'
% % % % % utemp1=utemp1(1:10)';
% % % % % utemp2=utemp2(1:10)';
% % % % 
% % % % % % uprevious = u(1:10,d-1)'
% % % % % % Qmat = Q(1:10)'
% % % % % % utemp1mat = BB*u(:,d-1);
% % % % % % utemp1=utemp1mat(1:10)'
% % % % % % utemp2mat = utemp2(1:10)'
% % % % 
% % % % % AA = full(AA);
% % % % 
% % % % 
% % % % % save mat A B Adebug Bdebug sizeM
% % % % % max(max(abs(AA - AAdebug)))
% % % % % max(max(abs(utemp2 - utemp2debug)))
% % % % % 
% % % % % a = max(max(abs(AA(1:sizeM,1:sizeM) - AAdebug(1:sizeM,1:sizeM)))) %a
% % % % % b = max(max(abs(AA(1:sizeM,sizeM+1:end) - AAdebug(1:sizeM,sizeM+1:end))))%b
% % % % % c = max(max(abs(AA(sizeM+1:end,1:sizeM) - AAdebug(sizeM+1:end,1:sizeM)))) %c
% % % % % d = max(max(abs(AA(sizeM+1:end,sizeM+1:end) - AAdebug(sizeM+1:end,sizeM+1:end)))) %d
% % % % 
% % % % % 
% % % % 
% % % % size(Q)
% % % % size(AA)
% % % % size(BB)
% % % % size(u(:,d-1))
% % % % 
% % % % u(:,d) = AA\(BB*u(:,d-1) + Q);
% % % % 
% % % % % udebug = AAdebug\utemp2debug; %from matrices in C
% % % % 
% % % % % uMATLAB  = u(1:10,d)'
% % % % % uC = udebug(1:10)'
% % % % 
% % % % % max(max(abs(u(:,d) - udebug)))
% % % % % 
% % % % % a = max(abs(u(1:sizeM,d) - udebug(1:sizeM))) %a
% % % % % b = max(abs(u(sizeM+1:end,d) - udebug(sizeM+1:end))) %a%b
% % % % % 
% % % % % u(sizeM+1,d)
% % % % % udebug(sizeM+1)
% % % % 
% % % % 
% % % % % error('er')
% % % % 
% % % % % 
% % % % % ufinal=u(1:10,d)'
% % % % % 
% % % % % error('debug');
% % % % uNEW = u(:,d);

end