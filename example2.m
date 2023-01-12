%Sylvester equation convection-diffusion
n = 1024;
k =26; % Numero di poli

%bs = 10;

%[A,B,U,V]=Palitta_example1(n); %A non simmetrica, B simmetrica
[A,B,U,V]=Palitta_example2(n);%A e B non simmetriche

 %rng(1)
 %U=[U,randn(n,bs)];
 %V=[V,randn(n,bs)];

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
dlmwrite('example2_det.dat',resdet,'\t');

options.poles="det2";
[Xu, Xv, resdet2] = Swapped_Sylvester_adaptive(A, B, U, V, options);
dlmwrite('example2_det2.dat',resdet2,'\t');


%options.poles="norm";
%[Xu, Xv, resnorm] = Swapped_Sylvester_adaptive(A, B, U, V,options);




[LA,UA]=lu(A);
[LB,UB]=lu(-B');
[Xu,Xv,reskpik]=kpik_sylv(A,LA,UA,-B',LB,UB,-U,V,(k+2)/2,1e-16);
dlmwrite('example2_kpik.dat',[[2:2:k+2]',reskpik],'\t');

%[Xu, Xv,resext] = Swapped_Sylvester(A,B,U,V,kron(ones(k/2),[0,inf]));
%resext=resext./norm(U*V','fro')




semilogy(max(resdet(:,1),resdet(:,2)),resdet(:,3),'b-');
hold on
semilogy(max(resdet2(:,1),resdet2(:,2)),resdet2(:,3), 'r-');
hold on
semilogy([2:2:k+2]',reskpik, 'k-');
hold on
%semilogy(resext, 'k+');
%hold on
%semilogy(resnorm,'g-');
legend('Drusk-Simon', 'nostro','kpik','ext','norm')

