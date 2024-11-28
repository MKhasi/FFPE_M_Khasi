function [A B C]=JacobiPol(xx,N,scale1,scale2,alpha,beta)
%Computing the Jacobi Polynomials

a(1)=(3+alpha+beta)*(4+alpha+beta)/(4*(2+alpha+beta));
b(1)=(beta^2-alpha^2)*(3+alpha+beta)/(4*(2+alpha+beta)*(2+alpha+beta));
c(1)=(alpha+1)*(beta+1)*(4+alpha+beta)/(2*(2+alpha+beta)*(2+alpha+beta));

%Jacobi Polynomials
A(1,:)=1+0.0.*xx;
A(2,:)=0.5*(alpha+beta+2)*(scale1*xx+scale2)+0.5*(alpha-beta);
B(1,:)=0.0.*xx;
B(2,:)=0.5*(alpha+beta+2)*scale1+0.0.*xx;
C(1,:)=0.0.*xx;
C(2,:)=0.0.*xx;

for n=2:N-1
    A(n+1,:)=(a(n-1)*(scale1*xx+scale2)-b(n-1)).*A(n,:)-c(n-1).*A(n-1,:);
    B(n+1,:)=scale1*(a(n-1).*A(n,:)+(a(n-1)*(scale1*xx+scale2)-b(n-1)).*B(n,:)-c(n-1).*B(n-1,:));
    C(n+1,:)=scale1^2*(2*a(n-1).*B(n,:)+(a(n-1)*(scale1*xx+scale2)-b(n-1)).*C(n,:)-c(n-1).*C(n-1,:));
    a(n)=(2*n+1+alpha+beta)*(2*n+2+alpha+beta)/(2*(n+1)*(n+1+alpha+beta));
    b(n)=(beta^2-alpha^2)*(2*n+1+alpha+beta)/(2*(n+1)*(n+1+alpha+beta)*(2*n+alpha+beta));
    c(n)=(alpha+n)*(beta+n)*(2*n+2+alpha+beta)/((n+1)*(n+1+alpha+beta)*(2*n+alpha+beta));
end