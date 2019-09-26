function [A BB] = coeff(yn,maturity,lamda,mu_q,omega)

%define the NS factor loading
c1 = ones(120, 1);

c2 = zeros(120, 1);

c3 = zeros(120, 1);

for i = 1:120
c1(i) = -i;
    
c2(i) = -(1 - exp(-lamda*i))/(lamda);
    
c3(i) = i*exp(-lamda*i)+c2(i);

end;

BB = [c1 c2 c3];

%define adjust term
A=zeros(120,1);
for i = 2:120
A(i,1)=A(i-1,1)+BB(i-1,:)*mu_q+1/2*BB(i-1,:)*omega*BB(i-1,:)';
end

%define regression coefficient
for i = 1:120
BB(i,1) = -BB(i,1)/i;
    
BB(i,2) = -BB(i,2)/i;
    
BB(i,3) = -BB(i,3)/i;

A(i,1)= -A(i,1)/i;

end;

A=A(maturity,:);
BB=BB(maturity,:);