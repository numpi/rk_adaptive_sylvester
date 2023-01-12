%Sylvester equation convection-diffusion
n = 1024;
k =26; % Maxit

[A,B,U,V]=matrices_example2(n);

options=[];
options.maxit=k;
options.tol=1e-16;
options.real=true;
options.mA=eigs(A,1,'smallestreal','Maxiterations',1e5);
options.MA=eigs(A,1,'largestreal','Maxiterations',1e5);
options.mB=eigs(B,1,'smallestreal','Maxiterations',1e5);
options.MB=eigs(B,1,'largestreal','Maxiterations',1e5);

options.poles="det";
[Xu, Xv, resdet] = rk_adaptive_sylvester(A, B, U, V,options);

%Uncomment to save data
%dlmwrite('example2_det.dat',resdet,'\t');

options.poles="det2";
[Xu, Xv, resdet2] = rk_adaptive_sylvester(A, B, U, V, options);

%Uncomment to save data
%dlmwrite('example2_det2.dat',resdet2,'\t');

semilogy(max(resdet(:,1),resdet(:,2)),resdet(:,3),'b-');
hold on
semilogy(max(resdet2(:,1),resdet2(:,2)),resdet2(:,3), 'r-');

legend('Drusk-Simon', 'det2')

function [A,B,U,V]=matrices_example2(n)
h=1/(n+1);
vi=0.0083;
w1=@(x)(1+(1+x).^2/4);
w2=@(y) (1/2)*y;

t=linspace(0,1,n);
Phi1=diag(w1(t));
Psi1=eye(n);
Phi2=eye(n);
Psi2=diag(w2(t));

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
T1=2*eye(n)-diag(ones(1,n-1),1)-diag(ones(1,n-1),-1);
T1=(1/h^2)*T1;
T2=T1';
B1=diag(ones(1,n-1),1)-diag(ones(1,n-1),-1);
B1=(1/(2*h))*B1;
B2=B1';
A=-vi*T1+Phi1*B1;
B=+vi*T2-B2*Psi2;
end


