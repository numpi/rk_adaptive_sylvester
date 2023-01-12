function [A,B,U,V]=Palitta_example2(n)
h=1/(n+1);
vi=0.0083;
%example1
w1=@(x)(1+(1+x).^2/4);
w2=@(y) (1/2)*y;

t=linspace(0,1,n);
Phi1=diag(w1(t));
Psi1=eye(n);
Phi2=eye(n);
Psi2=diag(w2(t));

[X,Y]=meshgrid(t);
%F=sinh((X+Y).^2)+cos((X-Y).^2);
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
%T1(1,1:2)=[1,0];
%T1(n,n-1:n)=[0,1];
T1=(1/h^2)*T1;
T2=T1';
B1=diag(ones(1,n-1),1)-diag(ones(1,n-1),-1);
%B1(1,2)=0;
%B1(n,n-1)=0;
B1=(1/(2*h))*B1;
B2=B1';
%boundary conditions
% y0=0*ones(n,1);
% x0=0*ones(1,n);
% y1=0*ones(n,1);
% x1=0*ones(1,n);
% 
% F1=vi*(T1+(1/h^2)*eye(n)+Psi1(1,1)*Phi1*B1)*y0;
% F2=vi*x0*(T2+(1/h^2)*eye(n)+Phi2(1,1)*B2*Psi2);
% F3=vi*(T1+(1/h^2)*eye(n)+Psi2(n,n)*Phi2*B2)*y1;
% F4=vi*x1*(T2+(1/h^2)*eye(n)+Phi1(n,n)*B1*Psi1);
% 
% 
% if norm(F1)~=0
%     U=[U,F1];
%     V=[V,eye(n,1)];
% end
% 
% if norm(F2)~=0
%     U=[U,eye(n,1)];
%     V=[V,F2'];
% end
% 
% if norm(F3)~=0
%     U=[U,F3];
%     V=[V,[zeros(n-1,1);1]];
% end
% 
% if norm(F4)~=0
%     U=[U,[zeros(n-1,1);1]];
%     V=[V,F4'];
% end

A=-vi*T1+Phi1*B1;
B=+vi*T2-B2*Psi2;

