%Lyapunov equation diffusion
n = 1024;
k =19; % Numero di poli

% bs = 10;
% 
% rng(1)
% U = randn(n,bs);
% V=-U;
% 
% A = spdiags(ones(n,1) * [-1 2 -1], -1:1, n, n);
% B=-A;
[A,B,U,V]=Diffusion(n);
V=-V;

options=[];
options.maxit=k;
options.tol=1e-16;
options.real=true;
options.mA=eigs(A,1,'smallestreal','Maxiterations',1e5);
options.MA=eigs(A,1,'largestreal','Maxiterations',1e5);
options.mB=eigs(B,1,'smallestreal','Maxiterations',1e5);
options.MB=eigs(B,1,'largestreal','Maxiterations',1e5);

options.poles="det";
[Xu, Xv, resdet] = Swapped_Sylvester_adaptive(A, B, U, V,options);
dlmwrite('example1_det.dat',resdet,'\t');
options.poles="det2";
[Xu, Xv, resdet2] = Swapped_Sylvester_adaptive(A, B, U, V, options);
dlmwrite('example1_det2.dat',resdet2,'\t');
%options.poles="norm";
%[Xu, Xv, resnorm] = Swapped_Sylvester_adaptive2(A, B, U, V, k,options);
%resnorm=resnorm./norm(U*V','fro');

[LA,UA]=lu(A);
[LB,UB]=lu(-B);
[Xu,Xv,reskpik]=kpik_sylv(A,LA,UA,-B',LB,UB,-U,V,(k+1)/2,1e-16);
dlmwrite('example1_kpik.dat',[[2:2:k+1]',reskpik],'\t');




params.m=k;
params.tol=1e-16;     
params.smin=min(eig(-A)); params.smax=max(eig(-A));
params.ch=0;       
params.period=1;   
[Z,resrk]=RKSM_Lyap_real(A,U,params);
resrk=resrk./norm(U*V','fro');

%[Xu, Xv,resext] = Swapped_Sylvester(A,B,U,V,kron(ones(1,k/2),[0,inf]));
%resext=resext./norm(U*V','fro');
%dlmwrite('example1_ext.dat',[[1:k]',resext'],'\t');

semilogy(max(resdet(:,1),resdet(:,2)),resdet(:,3),'b-');
hold on
semilogy(max(resdet2(:,1),resdet2(:,2)),resdet2(:,3), 'r-');
hold on
semilogy([2:2:k+1],reskpik, 'k-');
hold on
%semilogy(resext, 'k+');
%semilogy(resnorm,'g-');
%hold on
%semilogy(resrk, 'o-');
legend('Drusk-Simon', 'nostro','kpik','extended','resnorm', 'resrk')

