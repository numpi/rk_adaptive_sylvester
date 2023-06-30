%Lyapunov equation diffusion
n = 2^12; % Size of the problem
k = 201; % Maximum number of iterations

if ~exist('rat_krylov', 'file')
    error('RKTOOLBOX not found, please download it from ' + ... 
        'http://guettel.com/rktoolbox/');
end

[A,B,U,V] = Diffusion(n);
V=-V;
time=zeros(1,3);

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
%str="diffusion"+num2str(options.tol)+".txt";
%fileID = fopen(str,'w');
%fprintf(fileID,'&iter & residual & time (s) \n');
%fclose(fileID);

options.poles="det";
t = tic;
[Xu, Xv, resdet] = rk_adaptive_sylvester(A, B, U, V,options);
time(1) = toc(t);
% To save data of the residuals, uncomment the following line
%dlmwrite('example1_det.dat',resdet,'\t');
iter(1)=resdet(end,1);
residual(1) = norm(A*Xu*Xv'-Xu*Xv'*B-U*V', 'fro')/norm(U*V', 'fro');

% To save data of timings, uncomment the following lines
%fileID = fopen(str,'a');
%fprintf(fileID,' ADM & $%i$ & $%.2e$ & $%.2f$ \n',iter(1), residual(1) ,time(1));
%fclose(fileID);

options.poles="det2";
tic
[Xu, Xv, resdet2] = rk_adaptive_sylvester(A, B, U, V, options);
time(2)=toc;
% To save data for the residuals, un comment the following line
%dlmwrite('example1_det2.dat',resdet2,'\t');
iter(2)=resdet2(end,1);
residual(2) = norm(A*Xu*Xv'-Xu*Xv'*B-U*V', 'fro') / norm(U*V', 'fro');

% To save data of timings, uncomment the following lines
%fileID = fopen(str,'a');
%fprintf(fileID,' sADM & $%i$ & $%.2e$ & $%.2f$ \n',iter(2), residual(2) ,time(2));
%fclose(fileID);

t = tic;
options.poles="ext";
[Xu, Xv, resext] = rk_adaptive_sylvester(A, B, U, V, options);
time(3) = toc(t);
% To save data for the residuals, un comment the following line
%dlmwrite('example1_ext.dat',resext,'\t');
iter(3) = resext(end,1);
residual(3) = norm(A*Xu*Xv'-Xu*Xv'*B-U*V', 'fro')/norm(U*V', 'fro');

% To save data of timings, uncomment the following lines
%fileID = fopen(str,'a');
%fprintf(fileID,' ext & $%i$ & $%.2e$ & $%.2f$ \n',iter(3), residual(3) ,time(3));
%fclose(fileID);

semilogy(max(resdet(:,1),resdet(:,2)),resdet(:,3),'b-');
hold on
semilogy(max(resdet2(:,1),resdet2(:,2)),resdet2(:,3), 'r-');
hold on
semilogy(max(resext(:,1),resext(:,2)),resext(:,3),'k-');
legend('ADM', 'sADM','extended')


%---------------------------------

function [A,B,U,V]=Diffusion(n)
h=1/(n-1);
vi=1;

t=linspace(0,1,n);

[X,Y]=meshgrid(t);
F=1./((1+X+Y));
[U,S,V] = svd(F);
c=2;
while c<n
if S(c,c)<1e-10
S=S(1:c-1,1:c-1);
U=U(:,1:c-1);
V=V(:,1:c-1);
c=n;
end
c=c+1;
end
U=U*sqrt(S);
V=V*sqrt(S);

T1=spdiags(ones(n,1) * [-1 2 -1], -1:1, n, n);
T1=(1/h^2)*T1;
T2=T1';

A=-vi*T1;
B=vi*T2;
end

