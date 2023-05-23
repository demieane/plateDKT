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

AAmat = full(AA(1:10,1:10));
BBmat = full(BB(1:10,1:10));

% % % fileID = fopen('../c/OUTDATA_FEM_DEBUG.bin','rb')
% % % rowsUsol = fread(fileID,1,'int')
% % % colsUsol = fread(fileID,1,'int')
% % % 
% % % for i = 1:rowsUsol
% % %     for j = 1:colsUsol
% % %         AAdebug(i,j)=fread(fileID,1,'double');
% % % %         Adebug(i,j)=fread(fileID,1,'double');
% % %     end
% % % end
% % % 
% % % rowsUsol1 = fread(fileID,1,'int')
% % % colsUsol1 = fread(fileID,1,'int')
% % % 
% % % for i = 1:rowsUsol1
% % %     for j = 1:colsUsol1
% % %         utemp2debug(i,j)=fread(fileID,1,'double');
% % % %         Bdebug(i,j)=fread(fileID,1,'double');
% % %     end
% % % end

% save mat A B Adebug Bdebug sizeM
% max(max(abs(A - Adebug)))
% max(max(abs(B - Bdebug)))
% 
% a = max(max(abs(A(1:sizeM,1:sizeM) - Adebug(1:sizeM,1:sizeM)))) %a
% b = max(max(abs(A(1:sizeM,sizeM+1:end) - Adebug(1:sizeM,sizeM+1:end))))%b
% c = max(max(abs(A(sizeM+1:end,1:sizeM) - Adebug(sizeM+1:end,1:sizeM)))) %c
% d = max(max(abs(A(sizeM+1:end,sizeM+1:end) - Adebug(sizeM+1:end,sizeM+1:end)))) %d
% %%
% 
% a = max(max(abs(B(1:sizeM,1:sizeM) - Bdebug(1:sizeM,1:sizeM)))) %a
% b = max(max(abs(B(1:sizeM,sizeM+1:end) - Bdebug(1:sizeM,sizeM+1:end))))%b
% c = max(max(abs(B(sizeM+1:end,1:sizeM) - Bdebug(sizeM+1:end,1:sizeM)))) %c
% d = max(max(abs(B(sizeM+1:end,sizeM+1:end) - Bdebug(sizeM+1:end,sizeM+1:end)))) %d



Q = (1 - theta)*dt*G(:,d-1) + (theta)*dt*G(:,d);
% utemp2 = BB*u(:,d-1) + Q;


% % uprevious = u(1:10,d-1)'
% % Qmat = Q(1:10)'
% % utemp1mat = BB*u(:,d-1);
% % utemp1=utemp1mat(1:10)'
% % utemp2mat = utemp2(1:10)'

% AA = full(AA);


% save mat A B Adebug Bdebug sizeM
% max(max(abs(AA - AAdebug)))
% max(max(abs(utemp2 - utemp2debug)))
% 
% a = max(max(abs(AA(1:sizeM,1:sizeM) - AAdebug(1:sizeM,1:sizeM)))) %a
% b = max(max(abs(AA(1:sizeM,sizeM+1:end) - AAdebug(1:sizeM,sizeM+1:end))))%b
% c = max(max(abs(AA(sizeM+1:end,1:sizeM) - AAdebug(sizeM+1:end,1:sizeM)))) %c
% d = max(max(abs(AA(sizeM+1:end,sizeM+1:end) - AAdebug(sizeM+1:end,sizeM+1:end)))) %d

% error('er')
u(:,d) = AA\(BB*u(:,d-1) + Q);
% 
% u(1:10,d)'
% 
% error('debug');
uNEW = u(:,d);

end