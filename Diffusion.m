function [A,B,U,V]=Diffusion(n)
h=1/(n-1);
vi=1;
%example1

t=linspace(0,1,n);

[X,Y]=meshgrid(t);
%F=sin((X+Y).^3);
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


 %boundary conditions
%  y0=0*ones(n,1);
%  x0=0*ones(1,n);
%  y1=0*ones(n,1);
%  x1=0*ones(1,n);
%  
%  F1=-vi*(1/h^2)*y0;
%  F2=-vi*(1/h^2)*x0;
%  F3=-vi*(1/h^2)*y1;
%  F4=-vi*(1/h^2)*x1;
%  
%  
%  if norm(F1)~=0
%      U=[U,F1];
%      V=[V,eye(n,1)];
%  end
%  
%  if norm(F2)~=0
%      U=[U,eye(n,1)];
%      V=[V,F2'];
%  end
%  
%  if norm(F3)~=0
%      U=[U,F3];
%      V=[V,[zeros(n-1,1);1]];
%  end
%  
%  if norm(F4)~=0
%      U=[U,[zeros(n-1,1);1]];
%      V=[V,F4'];
%  end
A=-vi*T1;
B=vi*T2;

