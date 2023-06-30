%Sylvester equation convection-diffusion
n = 2^12;
k = 201; % Maxit

if ~exist('rat_krylov', 'file')
    error('RKTOOLBOX not found, please download it from ' + ... 
        'http://guettel.com/rktoolbox/');
end

[A,B,U,V]=matrices_example2(n);
time=zeros(1,3);
iter=zeros(1,3);
residual=zeros(1,3);

options=[];
options.maxit=k;
options.tol=1e-8;
options.real=true;
% Whenever estimates for the eigenvalues are available, they can be
% specified as options.
% options.mA=eigs(A,1,'smallestreal','Maxiterations',1e5);
% options.MA=eigs(A,1,'largestreal','Maxiterations',1e5);
% options.mB=eigs(B,1,'smallestreal','Maxiterations',1e5);
% options.MB=eigs(B,1,'largestreal','Maxiterations',1e5);

% To save data of timings, uncomment the following lines
%str="conv-diff"+num2str(options.tol)+".txt";
%fileID = fopen(str,'w');
%fprintf(fileID,'&iter & residual & time (s) \n');
%fclose(fileID);

options.poles="ADM";
tic
[Xu, Xv, resADM] = rk_adaptive_sylvester(A, B, U, V,options);
time(1)=toc;
iter(1)=resADM(end,1);
residual(1)=norm(A*Xu*Xv'-Xu*Xv'*B-U*V', 'fro')/norm(U*V', 'fro');

% To save data of timings, uncomment the following lines
%fileID = fopen(str,'a');
%fprintf(fileID,' ADM & $%i$ & $%.2e$ & $%.2f$ \n',iter(1), residual(1) ,time(1));
%fclose(fileID);

% To save data for the residuals, un comment the following line
%dlmwrite('example2_ADM.dat',resADM,'\t');



options.poles="sADM";
tic
[Xu, Xv, ressADM] = rk_adaptive_sylvester(A, B, U, V, options);
time(2)=toc;
iter(2)=ressADM(end,1);
residual(2)=norm(A*Xu*Xv'-Xu*Xv'*B-U*V', 'fro')/norm(U*V', 'fro');

% To save data of timings, uncomment the following lines
%fileID = fopen(str,'a');
%fprintf(fileID,' sADM & $%i$ & $%.2e$ & $%.2f$ \n',iter(2), residual(2) ,time(2));
%fclose(fileID);

% To save data for the residuals, un comment the following line
%dlmwrite('example2_sADM.dat',ressADM,'\t');

tic
options.poles="ext";
[Xu, Xv, resext] = rk_adaptive_sylvester(A, B, U, V, options);
time(3)=toc;
iter(3)=resext(end,1);
residual(3)=norm(A*Xu*Xv'-Xu*Xv'*B-U*V', 'fro')/norm(U*V', 'fro');

% To save data for the residuals, un comment the following line
%dlmwrite('example2_ext.dat',resext,'\t');

% To save data of timings, uncomment the following lines
%fileID = fopen(str,'a');
%fprintf(fileID,' ext & $%i$ & $%.2e$ & $%.2f$ \n',iter(3), residual(3) ,time(3));
%fclose(fileID);

semilogy(max(resADM(:,1),resADM(:,2)),resADM(:,3),'b-');
hold on
semilogy(max(ressADM(:,1),ressADM(:,2)),ressADM(:,3), 'r-');
hold on
semilogy(max(resext(:,1),resext(:,2)),resext(:,3),'k-');
legend('ADM', 'sADM','extended')


function [A,B,U,V]=matrices_example2(n)
h=1/(n+1);
vi=0.0083;
w1=@(x)(1+(1+x).^2/4);
w2=@(y) (1/2)*y;

t=linspace(0,1,n).';
Phi1=spdiags(w1(t), 0, n, n);
Psi2=spdiags(w2(t), 0, n, n);

[X,Y]=meshgrid(t);
F=1./((1+X+Y));
[U,S,V] = svd(F);
c=2;
while c>1
if S(c,c)<1e-10
S=S(1:c-1,1:c-1);
U=U(:,1:c-1);
V=V(:,1:c-1);
c=0;
end
c=c+1;
end
U=U*sqrt(S);
V=V*sqrt(S);
T1=spdiags(ones(n,1) * [-1 2 -1],-1:1, n, n);
T1=(1/h^2)*T1;
T2=T1';
B1=spdiags(ones(n,1) * [-1 1], [-1 1], n, n);
B1=(1/(2*h))*B1;
B2=B1';
A=-vi*T1+Phi1*B1;
B=+vi*T2-B2*Psi2;
end


