function [Z,resvals]=RKSM_Lyap_real(A,B,params)
%
%function Z=RKSM_Lyap_real(A,B,params)
%      
% Approximately Solve  
%                A X  +  X A'  + B B' = 0
%     X \approx Z Z'
%
% by the Rational Krylov subspace method 
% (Order reduction onto the Rational Krylov subspace)
%
% Input:

% A               coeff. matrix.  n x n
% B               rhs factor   n x p

% params.m        max space dimension allowed
% params.tol      Algebraic stopping tolerance 
% params.smin,params.smax estimates for real spectral interval
%                 associated with field of values of A'
%                 e.g., smax=norm(A,1); smin=smax/condest(A);
% params.ch       ch=1  complex poles  ch=0 real poles
% params.period   how often check convergence (period=1 means at each iteration)
%

% Hints:
% 2) Provide "comfortable" (loose bounds) estimates s1, emax

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% If this code is used, please cite:

% For the all-real version of the code:
%
% Kirsten, G., & Simoncini, V. (2019). Order reduction methods for solving 
% large-scale differential matrix Riccati equations. 
% arXiv preprint arXiv:1905.12119.
%
% For the overall rational Krylov Lyapunov eqn solver:
%
% V. Druskin & V.Simoncini 
% Adaptive rational Krylov subspaces for large-scale dynamical systems 
% Systems & Control Letters, 60 (2011), pp. 546-560. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


tic
   m=params.m;
   tol=params.tol;
   s1=params.smin;
   emax=params.smax;
   ch=params.ch;
   period=params.period;
   
   dim1 = [];

[n,n]=size(A);
B=full(B);
p=size(B,2);
I=speye(p);O=0*I;
In=speye(n);

csize = size(B,2);

[V,rr]=qr(B,0); 
nrmb=norm((rr),'fro')^2; 
beta=V'*B; beta2=beta*beta';
errtot=[];

VV=V;

H=sparse(p*(m+2),p*(m+1));
nrmrestotnew=[];
%nrma=norm(A,'fro');

if (norm(A-A',1)<1e-14), symm=1;
fprintf('The matrix A is symmetric \n')
else symm=0;
fprintf('The matrix A is nonsymmetric \n')
end

fprintf('\n')

fprintf('     # iter      rel.res.\n')

 newAv=A*V;
 %newAv=A'*V;
 K=full(V'*newAv); 
 s=s1(1);
 eH=eig(K);
 eHpoints = sort([s1(:)',emax]);
 snew=newpolei(eHpoints,eH,s1(1)*ones(p,1));
 if real(snew)<0, snew=-real(snew)+sqrt(-1)*imag(snew);end
 s=[s,snew];

% additional steps
cmplxflag=0;
itsinner=0;


i=0;

itcheck = i;
while i < m


  i=i+1;

  paired=0;
  itp=1;
  cflag = 0;
  Vwrk = V;

  while (paired==0),

    i1=i+1; it = 0; t=0.;
    w=Vwrk;
 
    
     its_cg=0; res_cg=0;
     wrk1 = (A-snew*In)\w; 
     
     %%%% All real basis implementation for RKSM %%%%%
     
     if imag(wrk1) ~= 0 & cflag == 0
         wrk = real(wrk1);
         cflag = 1;
     elseif imag(wrk1) ~= 0 & cflag == 1
         wrk = imag(wrk1);
     else
         wrk = wrk1;
     end
     
     
% Gram-Schmidt step
    jms=(i-1)*p+1;j1s=(i+1)*p;js=i*p;js1=js+1;
    Bm(jms:js,1:csize)= V'*B;
 
    for it=1:2,
      for kk=1:i
        k1=(kk-1)*p+1; k2=kk*p; 
        ww=VV(1:n,k1:k2);
        gamma=ww'*wrk;
        H(k1:k2,jms:js) = H(k1:k2,jms:js)+ gamma;
        wrk = wrk - ww*gamma;
      end
    end
    
    [V,hinv]=qr(wrk,0); H(js1:j1s,jms:js)=hinv; %hinv = inv(hinv);
    
    if (cmplxflag), 
        snew=conj(snew); s=[s,snew];cmplxflag=0;
        newAv=A*V;
        %newAv=A'*V;
        D = kron(spdiags(s(2:end)),I);
        g = VV'*newAv;
        g1 = g; 
        g2 = V'*A*VV; g3 = V'*A*V;
        %g2 = V'*A'*VV; g3 = V'*A'*V;
        K = [K g1; g2, g3];
        VV=[VV,V];
        i=i+1; itp=itp+1;
    else, 
        paired=1; 
    end
  end

  VVnew = [VV, V];
  ih1=i1; ih=i;
  newAv=A*V;
    %newAv=A'*V;
  D = kron(spdiags(s(2:end)),I);
  g = VV'*newAv;
    
    
  if (symm), K=(K+K')/2; end

  if (rem(i,period)==0)

% Solve the projected problem
      Y=lyap(K,Bm*Bm'); 

% computed residual   (exact, in exact arithmetic) cheaper computation possible
     u1=newAv-VV*g;
     d=-VV*(Y*(H(1:ih*p,1:ih*p)'\[sparse(p*(ih-1),p);I])*H(p*ih+1:p*ih1,p*ih-p+1:p*ih)');
     U=full([-V*s(end),  d u1 ]);
     rr=qr(U,0); rr=triu(rr(1:size(rr,2),:));
     

%backward error
     nrmres = norm(rr*sparse([O I O; I O I; O I O ])*rr','fro');

% relative residual norm
     nrmresnew = (nrmres)/nrmb;
 
     nrmrestotnew = [nrmrestotnew, nrmresnew];

     dim = size(VV,2);
     dim1 = [dim1,dim];

     disp([i,nrmresnew])

     if (nrmresnew<tol), 
        break
     end
  end 

% New poles and zeros
  eH=sort(eig(K));
  eHorig=eH;
  %eK=sort(eig(full(K)));

  if (ch)                     % Complex poles. Compute set for next complex pole of r_m

    if (any(imag(eH)) ~=0 & max(abs(imag(eH)))>1e-5 & length(eH)>2) % Roots lambdas come from convex hull too
     eH=[eH;-emax];
      ij=convhull(real(eH),imag(eH)); eH=eH(ij);
      ieH=length(eH); missing=ih*p-ieH;
      while missing>0,                         % include enough points from the border
        neweH=(eH(1:ieH-1)+eH(2:ieH))/2;missing=ih*p-length(eH);
        eH=[eH;neweH];
      end
    % eH=eH(1:ih);
      eHpoints=-eH;
      eH=eHorig;
    else                                  % if all real eigs, no convex hull possible
      eHpoints = sort([s1; emax.';-real(eH)]);
    end


  else   % Real poles s from real set. Compute complex roots of r_m via Ritz convex hull
     
     if (any(imag(eH)) ~=0 & length(eH)>2)    % Roots lambdas come from convex hull too
       eH=[eH;-s1;-emax.'];
       ij=convhull(real(eH),imag(eH)); eH=eH(ij);
       ieH=length(eH); missing=ih*p-ieH;
       while missing>0, % include enough points from the border
         neweH=(eH(1:ieH-1)+eH(2:ieH))/2;
         eH=[eH;neweH];
         missing=ih*p-length(eH);
       end
       eH=eH(1:ih*p);
     end
      eHpoints = sort([s1; emax.';-real(eH)]);
      eH=eHorig;
  end


  gs=kron(s(2:end),ones(1,p))';

  snew = newpolei(eHpoints,eH,gs);
  if real(snew)<0, snew=-real(snew)+sqrt(-1)*imag(snew);end


% If pole is complex, include its conjugate

  if (imag(snew) ~=0), cmplxflag=1;end
  s=[s,snew];

  g1 = g; 
  g2 = V'*A*VV; g3 = V'*A*V;
  %g2 = V'*A'*VV; g3 = V'*A'*V;
  K = [K g1; g2, g3];
  VV=[VV,V];
    
    
end


% factored solution 

[uY,sY]=eig(Y); [sY,id]=sort(diag(sY));
sY=flipud(sY); uY=uY(:,id(end:-1:1));
is=sum(abs(sY)>1e-8);
Y0 = uY(:,1:is)*diag(sqrt(sY(1:is))); 
Z = VV(:,1:size(Y0,1))*Y0; 
final_rank=is;
t2=toc;
 fprintf(' \n')
fprintf('Space dim %d  Solution rank %d time %d\n',size(VV,2),is,t2);
 
resvals=nrmrestotnew*nrmb;
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r=ratfun(x,eH,s)
 
for j=1:length(x)
r(j)=abs(prod( (x(j)-s)./(x(j)-eH) ));
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r=ratfuni(x,eH,s)

for j=1:length(x)
r(j)=1./abs(prod( (x(j)-s)./(x(j)-eH) ));
end

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function snew=newpolei(eHpoints,eH,s)

for j=1:length(eHpoints)-1
   %snew(j) = fminbnd( @(x) ratfuni(x,eH,s), eHpoints(j),eHpoints(j+1),optimset('TolX',1e-3));
    sval=linspace(eHpoints(j),eHpoints(j+1),20);
   %sval=linspace(eHpoints(j),eHpoints(j+1),100);
    [sf,jx] = max (abs(ratfun(sval,eH,s)));
    snew(j)=sval(jx);
end
[sn,jx]=max(abs(ratfun(snew,eH,s)));
snew=snew(jx);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function RA=compute_RA(A,eH,s);
%
%ne=length(eH);
%s=[1e20,s];
%I=speye(size(A));
%RA=I;
%for k=1:ne,
%   RA = (A-eH(k)*I)/(A-s(k)*I)*RA;
%end
%return