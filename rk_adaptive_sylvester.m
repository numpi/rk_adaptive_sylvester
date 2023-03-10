function [Xu, Xv,resval] = rk_adaptive_sylvester(A,B,u,v, options)
%[Xu,Xv,resval] = rk_adaptive_sylvester(A,B,u,v,options)
%
% Rational Krylov method for Sylvester equations that ensures the last pole
% is equal to infinity at each step, to check the residual. 
%
% The function solves the Sylvester equation
% 
%    A X - X B - u v' = 0
%
% assuming sep(A,B) >0
%
%Input parameters:
%  A, B, square matrices of dimension nA x nA and nB x nB, respectively;
%  u,v block vectors of dimension nA x bs and nB x bs, respectively;
% 
%  options.maxit = max space dimension, defalut min(sqrt(size(A)),sqrt(size(B)));
%
%  options.tol = max final accuracy (in terms of relative residual),
%  default 1e-10.
%
%  Stopping criterion:
%       ||A Xu Xv' - Xu Xv' B - u v'||_F
%       --------------------------------------  < tol
%                 || u  v'||_F
%
%  options.poles determines how to adaptively chose poles. The possible
%  choices are options.poles='det', options.poles='det2', default is 'det2'
%
%  The pole selection algorithm are described in detail in [1]. Algorithm
%  'det' is the heuristic proposed in [2].
%
%
%  options.mA, option.MA determines a lower and upper bound, respectively, for
%  the real part of the eigenvalues of A;
%
%  options.mB, option.MB determines a lower and upper bound, respectively, for
%  the real part of thr eigenvalues of B;
%
%  options.real='true' runs the algorithm for A and B real;
%
% Output parameters:
%
%  Xu,Xv = solution factors   X_approx = Xu Xv'
%  resval   history of the vector
%   [number of iteration A, number of iterations B, relative residual norm]
%
% References:
% [1] Casulli, A, Robol, L., An efficient block rational Krylov solver for Sylvester
% equations with adaptive pole selection.
% [2] Druskin V., Simoncini V.,Adaptive rational krylov subspaces for
% large-scale dynamical systems. Systems & Control Letters,2011.


if nargout==3
    resval=[];
end

if nargin==4
    options=[];
end

if isfield(options, 'tol')==0
    options.tol=1e-10;
end

if isfield(options, 'maxit')==0
    options.maxit=min(sqrt(size(A,1)),sqrt(size(B,1)));
end
maxit=options.maxit;

if isfield(options, 'poles')==0
    options.poles="det2";
end

normUV=norm(u*v','fro');

[VA,KA,HA] = rat_krylov(A, u, inf);
bs=size(u,2);
Ap=HA(1:bs,1:bs)/KA(1:bs,1:bs);

[VB,KB,HB] = rat_krylov(B', v, inf);
Bp=HB(1:bs,1:bs)/KB(1:bs,1:bs);

C=(VA(:,1:bs)'*u)*(VB(:,1:bs)'*v)';

if isfield(options, 'mA')==0
    options.mA=eigs(A,1,'smallestreal');
end
mA=options.mA;
if isfield(options, 'MA')==0
    options.MA=eigs(A,1,'largestreal');
end
MA=options.MA;

if isfield(options, 'mB')==0
    options.mB=eigs(B,1,'smallestreal');
end
mB=options.mB;
if isfield(options, 'MB')==0
    options.MB=eigs(B,1,'largestreal');
end
MB=options.MB;

if abs(MB)<abs(mB)
    snewA=MB;
else
    snewA=mB;
end

if ~isreal(snewA)
snewA=[snewA,conj(snewA)];
end

if abs(MA)<abs(mA)
    snewB=MA;
else
    snewB=mA;
end

if ~isreal(snewB)
snewB=[snewB,conj(snewB)];
end

sA=[MB+mB-snewA(1),snewA];
sB=[MA+mA-snewB(1),snewB];
k=1;
h=1;
while k<=maxit && h<=maxit
        k=k+1;
        h=h+1;
    if options.real
        [VA,KA,HA,Ap] = Swapped_update_real(A,VA,KA,HA,snewA, Ap);
        [VB,KB,HB,Bp] = Swapped_update_real(B',VB,KB,HB,snewB, Bp);
        C(end+length(snewA)*bs,end+length(snewB)*bs)=0;
        k=k+length(snewA)-1;
        h=h+length(snewB)-1;
    else
        [VA,KA,HA,Ap] = Swapped_update(A,VA,KA,HA,snewA, Ap);
        [VB,KB,HB,Bp] = Swapped_update(B',VB,KB,HB,snewB, Bp);
        C(end+bs,end+bs)=0;
    end
    
    
    [QA,TA]=schur(Ap, 'complex');
    [QB,TB]=schur(Bp, 'complex');
    
    Xp = matlab.internal.math.sylvester_tri(TA, -TB, QA' * C * QB, 'I', 'I', 'transp');
    
    Y = QA * Xp * QB';
    
    if options.poles=="det"
        snewA = newpole_adaptive_det(mB,MB,diag(TB),diag(TA),sA);
        snewB = newpole_adaptive_det(mA,MA,diag(TA),diag(TB),sB);
    elseif options.poles=="det2"
        snewA = newpole_adaptive_det2(mB,MB,diag(TB),diag(TA),sA);
        snewB = newpole_adaptive_det2(mA,MA,diag(TA),diag(TB),sB);
    end
    
    if options.real
        if ~isreal(snewA)
            snewA(2)=conj(snewA);
        end
        if ~isreal(snewB)
            snewB(2)=conj(snewB);
        end
    end
    
    sA = [sA, snewA];
    sB = [sB, snewB];
    
    
    
    resA = norm( ( HA(end-bs+1:end, :) / KA(1:end-bs, :) ) * Y , 'fro');
    resB = norm( ( HB(end-bs+1:end, :) / KB(1:end-bs, :) ) * Y' , 'fro');
    
    % 2-norm
    %resnrm = max(resA, resB);
    %fprintf('Iteration %d, res = %e\n', k, resnrm);
    
    %Frobenius norm
    resnrm = hypot(resA, resB)/normUV;
    fprintf('IterationA %d,IterationB %d, res = %e\n',k,h, resnrm);
    if resnrm<options.tol
        break
    end
    if nargout==3
        resval=[resval;[k,h,resnrm]];
    end
end
if k>= maxit
    fprintf('Stopped at iteration %d\n', k);
end
Xu = VA(:,1:end-bs) * Y;
Xv = VB(:,1:end-bs);
end


function [V,K,H,Ap] = Swapped_update(A,V,K,H,s, Ap)
%Step of rational Krylov method to add the pole s
% and swap poles ensuring that that the last pole
% is equal to infinity.

bs = size(H, 1) - size(H, 2);
param.extend=bs;

    [V,K,H] = rat_krylov(A,V,K,H,s,param);
    k=size(K,2);
    [G,K(end-2*bs+1:end,end-bs+1:end)]=qr(K(end-2*bs+1:end,end-bs+1:end));
    G=G';
    V(:,end-2*bs+1:end)=V(:,end-2*bs+1:end)*G';
    H(end-2*bs+1:end,end-2*bs+1:end)=G*H(end-2*bs+1:end,end-2*bs+1:end);
    
    [G,~] = qr(H(end-bs+1:end, end-2*bs+1:end).');
    G=G(:,end:-1:1);
    H(:, end-2*bs+1:end)=H(:, end-2*bs+1:end)*G;
    K(:, end-2*bs+1:end)=K(:, end-2*bs+1:end)*G;
      
    e=zeros(k,bs);
    e(end-bs+1:end,end-bs+1:end)=eye(bs);
    Ap(1:end+bs,end+1:end+bs)=H(1:end-bs,1:end)*(K(1:end-bs,1:end)\e);
    Ap(end-bs+1:end,1:end)=H(end-2*bs+1:end-bs,1:end)/K(1:end-bs,1:end);
 
    
    
end


function [V,K,H,Ap] = Swapped_update_real(A,V,K,H,s, Ap)
%Step of rational Krylov method to add the pole s (and its conjugate if s
%is not real), and swap poles ensuring that that the last pole
% is equal to infinity.

bs = size(H, 1) - size(H, 2);
param.extend=bs;
param.real=1;

    [V,K,H] = rat_krylov(A,V,K,H,s,param);
    k=size(K,2);
    
    for i=length(s):-1:1
    [G,~]...
        =qr(K(end-(i-1)*bs-2*bs+1:end-(i-1)*bs,end-(i-1)*bs-bs+1:end-(i-1)*bs));
    G=G';
    V(:,end-(i-1)*bs-2*bs+1:end-(i-1)*bs)=V(:,end-(i-1)*bs-2*bs+1:end-(i-1)*bs)*G';
    K(end-(i-1)*bs-2*bs+1:end-(i-1)*bs,end-(i-1)*bs-2*bs+1:end)...
        =G*K(end-(i-1)*bs-2*bs+1:end-(i-1)*bs,end-(i-1)*bs-2*bs+1:end);
    H(end-(i-1)*bs-2*bs+1:end-(i-1)*bs,end-(i-1)*bs-2*bs+1:end)...
        =G*H(end-(i-1)*bs-2*bs+1:end-(i-1)*bs,end-(i-1)*bs-2*bs+1:end);
    
    [G,~] = qr(H(end-(i-1)*bs-bs+1:end-(i-1)*bs, end-(i-1)*bs-2*bs+1:end-(i-1)*bs).');
    G=G(:,end:-1:1);
    H(:, end-(i-1)*bs-2*bs+1:end-(i-1)*bs)=H(:, end-(i-1)*bs-2*bs+1:end-(i-1)*bs)*G;
    K(:, end-(i-1)*bs-2*bs+1:end-(i-1)*bs)=K(:, end-(i-1)*bs-2*bs+1:end-(i-1)*bs)*G;

    e=zeros(k-(i-1)*bs,bs);
    e(end-bs+1:end,:)=eye(bs);
    Ap(1:end+bs,end+1:end+bs)=...
        H(1:end-(i-1)*bs-bs,1:end-(i-1)*bs)*(K(1:end-(i-1)*bs-bs,1:end-(i-1)*bs)\e);
    Ap(end-bs+1:end,1:end)=...
        H(end-(i-1)*bs-2*bs+1:end-(i-1)*bs-bs,1:end-(i-1)*bs)/K(1:end-(i-1)*bs-bs,1:end-(i-1)*bs);
    end
end


function np = newpole_adaptive_det(a,b,eigenvaluesA,eigenvaluesB, poles)
%Computes the newpole for rational Krylov method maximizing the determinant. 

poles = poles(:);

eigenvaluesB = eigenvaluesB(:);

bs = length(eigenvaluesB) / length(poles);
poles = kron(ones(bs, 1), poles);

if isreal(eigenvaluesA)  
    eHpoints=sort([a;b;eigenvaluesA]);
    maxvals=zeros(2,length(eHpoints)-1);
    for j=1:length(eHpoints)-1
    t=linspace(eHpoints(j),eHpoints(j+1),20);
    [maxvals(1,j),jx]= max(abs(ratfun(t, eigenvaluesB, poles)));
    maxvals(2,j)=t(jx);
    end
    [~,jx]=max(maxvals(1,:));
    np=maxvals(2,jx);
else
    x=[eigenvaluesA;a;b];
    k = convhull(real(x),imag(x));
    maxvals=zeros(2,length(k)-1);
    for i=1:length(k)-1
        t=linspace(x(k(i)),x(k(i+1)),500);
        [maxvals(1,i),jx] =max(abs(ratfun(t, eigenvaluesB, poles)));
        maxvals(2,i)=t(jx);
    end
    [~,jx]=max(maxvals(1,:));
    np=maxvals(2,jx);
end
end


function np = newpole_adaptive_det2(a,b,eigenvaluesA, eigenvaluesB, poles)
%Computes the newpole for rational Krylov method maximizing the product
%of a selected set of eigenvalues. 

poles = poles(:);

eigenvaluesB = eigenvaluesB(:);
bs = length(eigenvaluesB) / length(poles);

if isreal(eigenvaluesA)
    eHpoints=sort([a;b;eigenvaluesA]);
    maxvals=zeros(2,length(eHpoints)-1);
    for i=1:length(eHpoints)-1
        t=linspace(eHpoints(i),eHpoints(i+1),20);
        vals=zeros(1,length(t));
        for j=1:length(t)
            [~,I]=sort(abs(t(j)-eigenvaluesB), 'ascend');
            sorteig=eigenvaluesB(I);
            vals(j)=abs(prod( (t(j)-poles(2:end))./(t(j)-sorteig(bs+1:bs:end)) )...
                /((t(j)-sorteig(1))));
        end
        [maxvals(1,i),jx]= max(vals);
        maxvals(2,i)=t(jx);
    end
    [~,jx]=max(maxvals(1,:));
    np=maxvals(2,jx);
    
else
    x=[eigenvaluesA;a;b];
    k = convhull(real(x),imag(x));
    maxvals=zeros(2,length(k)-1);
    for i=1:length(k)-1
        t=linspace(x(k(i)),x(k(i+1)),500);
        vals=zeros(1,length(t));
        for j=1:length(t)
            [~,I]=sort(abs(t(j)-eigenvaluesB), 'ascend');
            sorteig=eigenvaluesB(I);
            vals(j)=abs(prod( (t(j)-poles(2:end))./(t(j)-sorteig(bs+1:bs:end)) )...
                /((t(j)-sorteig(1))));
        end
        [maxvals(1,i),jx] = max(vals);
        maxvals(2,i)=t(jx);
    end
    [~,jx]=max(maxvals(1,:));
    np=maxvals(2,jx);
end

end

function r=ratfun(x,eH,s)

r=zeros(size(x));
 
for j=1:length(x)
r(j)=abs(prod( (x(j)-s(2:end))./(x(j)-eH(2:end)) )/((x(j)-eH(1))));
end

return

end