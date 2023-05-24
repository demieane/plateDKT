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

% sizeBB = size(BB);

Q = (1 - theta)*dt*G(:,d-1) + (theta)*dt*G(:,d);
Qtemp = Q(1:10)'
Qtemp1 = Q(sizeM-4:sizeM+10)'

u(:,d) = AA\(BB*u(:,d-1) + Q);
utemp2 = (BB*u(:,d-1) + Q);

utemp1  = BB*u(:,d-1);

utemp11 = utemp1(1:10)'
utemp22 = utemp2(1:10)'
% uprevious = u(1:10,d-1);

uNEW = u(:,d);





% % % % sizeM=size(Mglob,1);
% % % % 
% % % % A = [Mglob, sparse(sizeM,sizeM); sparse(sizeM,sizeM), speye(sizeM,sizeM)];
% % % % % B = -[C, Kglob; -speye(sizeM,sizeM), sparse(sizeM,sizeM)];
% % % % B = [C, Kglob; -speye(sizeM,sizeM), sparse(sizeM,sizeM)];
% % % % 
% % % % %lamda (1: implicit Euler, 1/2: Crank - Nicolson)
% % % % 
% % % % AA =  A + theta*dt*B;
% % % % BB =  A - (1 - theta)*dt*B;
% % % % 
% % % % % AA =  A - theta*dt*B;
% % % % % BB =  A + (1 - theta)*dt*B;
% % % % 
% % % % AAmat = full(AA(1:10,1:10));
% % % % BBmat = full(BB(1:10,1:10));
% % % % 
% fileID = fopen('../c/OUTDATA_FEM_DEBUG.bin','rb')
% rowsUsol = fread(fileID,1,'int')
% colsUsol = fread(fileID,1,'int')
% 
% for i = 1:rowsUsol
%     for j = 1:colsUsol
% %         udebug(i,j)=fread(fileID,1,'double');
% %         Qdebug(i,j)=fread(fileID,1,'double');
%         AAdebug(i,j)=fread(fileID,1,'double');
% %         Adebug(i,j)=fread(fileID,1,'double');
%     end
% end
% 
% rowsUsol1 = fread(fileID,1,'int')
% colsUsol1 = fread(fileID,1,'int')
% 
% for i = 1:rowsUsol1
%     for j = 1:colsUsol1
% %         utemp2debug(i,j)=fread(fileID,1,'double');
%         BBdebug(i,j)=fread(fileID,1,'double');
% %         Bdebug(i,j)=fread(fileID,1,'double');
%     end
% end
% % % % % 
% % % % % 
% % % % % 
% max(max(abs(AA - AAdebug)))
% max(max(abs(BB - BBdebug)))
% 
% fclose(fileID)
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