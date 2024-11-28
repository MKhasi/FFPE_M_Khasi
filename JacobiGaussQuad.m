function [X W]=JacobiGaussQuad(N,alpha,beta)
%computing the nodes and weights of the general Jacobi-Gauss quadrature
%base on "Shen, Jie, Tao Tang, and Li-Lian Wang. Spectral methods: algorithms, analysis and applications. Vol. 41. Springer Science & Business Media, 2011"
%Page 83

A=zeros(N,N);
A(1,1)=(beta-alpha)/(alpha+beta+2);
for i=2:N
    A(i,i)=(beta^2-alpha^2)/((2*i-2+alpha+beta)*(2*i+alpha+beta));
end
for i=2:N
    A(i-1,i)=sqrt(4*(i-1)*(i-1+alpha)*(i-1+beta)*(i-1+alpha+beta)/((2*i-3+alpha+beta)*(2*i-2+alpha+beta)^2*(2*i-1+alpha+beta)));
    A(i,i-1)=A(i-1,i);
end

[P D]=eig(A);
X=diag(D);
for i=1:N
    W(i,1)=(2^(alpha+beta+1)*gamma(alpha+1)*gamma(beta+1)*P(1,i)^2)/(gamma(alpha+beta+2));
end