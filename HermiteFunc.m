function [A B C]=HermiteFunc(xx,N,scale1,scale2,Flag)
%Computing Hermite Functions
%input: xx   an array
%output: A, B, C   a matrix which compute the Hermite Polynomials and there
%derivatives A_{ij}=H_i(x_j)
% if Flage=='D' compute the derivatives

u=scale1*xx+scale2;
A(1,:)=pi^(-0.25).*exp(-0.5*u.^2);
A(2,:)=pi^(-0.25).*exp(-0.5*u.^2).*sqrt(2).*u;
for i=3:N
  A(i,:)=sqrt(2/(i-1)).*u.*A(i-1,:)-sqrt((i-2)/(i-1)).*A(i-2,:);
end

if Flag=='D'
  %First derivative of Hermite Polynomials
  B(1,:)=-scale1.*u.*pi^(-0.25).*exp(-0.5*u.^2);
  B(2,:)=scale1.*pi^(-0.25).*exp(-0.5.*u.^2).*sqrt(2).*(1-u.^2);
  
  for i=3:N
    B(i,:)=scale1.*(sqrt(2/(i-1)).*(A(i-1,:)+u.*B(i-1,:))-sqrt((i-2)/(i-1)).*B(i-2,:));
  end
  
  %Second derivative of Hermite Polynomials
  C(1,:)=scale1.^2.*pi.^(-0.25).*exp(-0.5*u.^2).*(u.^2-1);
  C(2,:)=scale1.^2.*pi.^(-0.25).*exp(-0.5*u.^2).*sqrt(2).*(u.^3-3.*u);
  for i=3:N
    C(i,:)=scale1.^2.*(sqrt(2/(i-1)).*(2*B(i-1,:)+u.*C(i-1,:))-sqrt((i-2)/(i-1)).*C(i-2,:));
  end
end
