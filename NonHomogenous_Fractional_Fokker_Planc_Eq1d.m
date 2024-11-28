% This code is for Solving the following 1d Fractional Fokker–Planck Equations
% D^{a\lpha} u= a_1 u_xx+b_1 u_x+c_1 u+f
% u(x,0)=g(x)
% u(+-\infty,t)->h_1(t),h_2(t)
% by Hermite–Galerkin spectral method in space and Petrov–Galerkin method based on generalized
% Jacobi functions in time
% Written by Manoochehr Khasi - 11/13/2024

clc;
clear all;


EX='Ex1';

N=20;         %The number of Space grid points
M=2*N;        %The number of Time grid points
MM=zeros(N,N);
for i=1:N-1
  MM(i,i+1)=sqrt(0.5*i);
  MM(i+1,i)= MM(i,i+1);
end
x=eig(MM);

switch EX
  case 'Ex1'
    alpha=0.5;
    TM=1;
    a_1=1;
    b_1=0;
    c_1=0;
    MD=10;
    h_1=@(t)  t.^2;
    h_1alpha=@(t)   2.*t.^(2-alpha)./gamma(3-alpha);
    h_1alphat=@(t)  2.*t.^(1-alpha)./gamma(2-alpha);
    h_1t=@(t)  2.*t;
    h_2=@(t)  t.^2;
    h_2alpha=@(t)  2.*t.^(2-alpha)./gamma(3-alpha);
    h_2alphat=@(t) 2.*t.^(1-alpha)./gamma(2-alpha);
    h_2t=@(t)  2.*t;
    u=@(x,t) (t.^alpha+1).*exp(-x.^2)+t.^2;
    ut=@(x,t) (alpha.*t.^(alpha-1)).*exp(-x.^2)+2.*t;
    ux=@(x,t) -2.*x.*(t.^alpha+1).*exp(-x.^2);
    uxt=@(x,t) -2.*x.*(alpha.*t.^(alpha-1)).*exp(-x.^2);
    uxx=@(x,t) (4.*x.^2-2).*(t.^alpha+1).*exp(-x.^2);
    uxxt=@(x,t) (4.*x.^2-2).*(alpha.*t.^(alpha-1)).*exp(-x.^2);
    ualpha=@(x,t) (gamma(alpha+1)).*exp(-x.^2)+2.*t.^(2-alpha)./gamma(3-alpha);
    ualphat=@(x,t) 2.*t.^(1-alpha)./gamma(2-alpha)+0.*x.*t;
    W=@(x,t)  (1/(2*MD))*(h_2(t)-h_1(t)).*x+0.5*(h_2(t)+h_1(t));
    Wt=@(x,t)  (1/(2*MD))*(h_2t(t)-h_1t(t)).*x+0.5*(h_2t(t)+h_1t(t));
    Walpha=@(x,t)  (1/(2*MD))*(h_2alpha(t)-h_1alpha(t)).*x+0.5*(h_2alpha(t)+h_1alpha(t));
    Walphat=@(x,t)  (1/(2*MD))*(h_2alphat(t)-h_1alphat(t)).*x+0.5*(h_2alphat(t)+h_1alphat(t));
    g=@(x) u(x,0)-W(x,0);
    ff=@(x,t)  ualpha(x,t)-a_1.*uxx(x,t)-b_1.*ux(x,t)-c_1.*u(x,t);
    fft=@(x,t)  ualphat(x,t)-a_1.*uxxt(x,t)-b_1.*uxt(x,t)-c_1.*ut(x,t);
    f=@(x,t) ff(x,t)+c_1*W(x,t)+b_1.*(1/(2*MD))*(h_2(t)-h_1(t))-Walpha(x,t);
    ft=@(x,t) fft(x,t)+c_1*Wt(x,t)+b_1.*(1/(2*MD))*(h_2t(t)-h_1t(t))-Walphat(x,t);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Creating Matrix A
A=zeros(N,N);
for n=1:N
  A(n,n)=a_1*(-n+0.5)+c_1;
end

for n=1:N-1
  A(n,n+1)=b_1*sqrt(0.5*n);
end

for n=2:N
  A(n,n-1)=-b_1*sqrt(0.5*(n-1));
end

for n=1:N-2
  A(n,n+2)=0.5*a_1*sqrt(n*(n+1));
end

for n=3:N
  A(n,n-2)=0.5*a_1*sqrt((n-2)*(n-1));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
psim=zeros(N,1);
for k=1:floor(N/2)
  psim(2*k,1)=0.0;
  for l=1:k-1
    psim(2*k,1)=psim(2*k,1)+log((2*l+1)./(2*l));
  end
  psim(2*k,1)=sqrt(exp(psim(2*k,1)))*2*pi^0.25;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[PHI PHIx PHIxx]=HermiteFunc(x',N,1,0,'D');
[P D]=eig(A);
t=linspace(0,TM,M);

[TJ WJ]=JacobiGaussQuad(M,0,alpha);
TJ=0.5*TM+0.5*TM*TJ;
[PPol PPol1 PPol2]=JacobiPol(TJ',M,2/TM,-1,-alpha,alpha);
[PPolt PPolt1 PPolt2]=JacobiPol(t,M,2/TM,-1,-alpha,alpha);
[TL WL]=JacobiGaussQuad(M,0,0);     %for Gauss-Legendre Quadrature
TL=0.5*TM+0.5*TM*TL;
[LPol LPol1 LPol2]=JacobiPol(TL',M,2/TM,-1,0,0);
[LPolTJ LPolTJ1 LPolTJ2]=JacobiPol(TJ',M,2/TM,-1,0,0);
G=(PHI*diag(1./(PHI(N,:)'.^2))*g(x))/N;

xx=repmat(x,1,M);         TLTL=repmat(TL',N,1);
F=(PHI*diag(1./(PHI(N,:)'.^2)))*f(xx,TLTL)/N;
Gbar=P\G;
F=P\F;

fm=zeros(M,1);
fmp=zeros(M,1);
fmp(1)=1/(alpha*gamma(alpha));
fm(1)=0.5*fmp(1);
for j=2:M
  fmp(j)=fmp(j-1)*((j-1)/(j-1+alpha));
  fm(j)=(j-0.5)*fmp(j);
end

FF=(F*diag(WL)*LPol')*diag(fm);
FF(:,1)=FF(:,1)+(diag(D).*Gbar)/gamma(1+alpha);
Bq=(0.5*TM)^alpha*(PPol*diag(WJ)*LPolTJ')*diag(fm);
[Q TMat]=schur(Bq','complex');

dd=zeros(N,M);
FF=FF*Q;
for i=1:N
  dd(i,:)=FF(i,:)/(eye(M,M)-D(i,i)*TMat');
end
MU=dd*Q'*PPolt*diag(t.^alpha)+repmat(Gbar,1,M);
nummeshX=1000;    nummeshT=100*TM;
X=linspace(-MD,MD,nummeshX)';
T=linspace(0,TM,nummeshT);
[PPolT PPolT1 PPolT2]=JacobiPol(T,M,2/TM,-1,-alpha,alpha);
MUT=dd*Q'*PPolT*diag(T.^alpha)+repmat(Gbar,1,nummeshT);

[xmesh tmesh]=meshgrid(x,t);
[Xmesh Tmesh]=meshgrid(X,T);
PHIX=HermiteFunc(X',N,1,0,'NoD');
ZZ=PHIX'*P*MUT+W(Xmesh',Tmesh');
semilogy(X,abs(ZZ(:,nummeshT)-u(X,TM)));
hold on;
%Global Plots
SURF(Xmesh,Tmesh,real(u(Xmesh,Tmesh)),'x','t','u(x,t)');             %plot of Exact solution
SURF(Xmesh,Tmesh,real(ZZ)','x','t','U(x,t)');                        %plot of Approximate solution
SURF(Xmesh,Tmesh,abs(real(ZZ'-u(Xmesh,Tmesh))),'x','t','Error');     %plot of the Error
title(['ERR=',num2str(max(max(abs(real(ZZ'-u(Xmesh,Tmesh))))))]);   
fprintf('\n   N= %d  %2.4e  %2.4e; ', N,max(abs(ZZ(:,nummeshT)-u(X,TM))),norm(ZZ(:,nummeshT)-u(X,TM)));
