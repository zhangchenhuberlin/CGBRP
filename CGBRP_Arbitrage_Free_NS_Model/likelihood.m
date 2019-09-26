function [s2] = likelihood(y,x,A,BB,h)
s2=0;
for i=1:size(y,1)
s2=s2+(y(i,:)'- BB*x(:,i) -A)'*(y(i,:)'- BB*x(:,i)-A);
end
s2=-s2*0.5*h;