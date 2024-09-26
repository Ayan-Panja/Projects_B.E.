clc
clear all
h=0.001;
z=0.1:h:2;
n=length(z);
R=zeros(1,n);
%R(1)=0;
%z(1)=0;
f=@(z) (2-z.^2)/((4*z.^2-z.^4).^(1/2));
for i=1:n-1
    k1=h*f(z(i));
    k2=h*f(z(i)+h/2);
    k4=h*f(z(i)+h);
    k=(k1+4*k2+k4)/6;
    R(i+1)=R(i)+k;
end

