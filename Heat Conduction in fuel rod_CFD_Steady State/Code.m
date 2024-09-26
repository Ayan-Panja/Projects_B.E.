close all;
clc;
%Specifying no of grids
n=201;
nx=200; 
ny=200;
Y=0.01;
%Grid size same for both
X=Y;
%Giving input value of Biot no.
Bi=input('enter Biot no:  '); 
%Initializing the variables
T=zeros(n); 
itr=0;
%Implementing Gauss-Siedel method
while itr<=500
       itr=itr+1;
        Told=T;
          %Insulated boundary conditions
             for j=102:200
                  T(101,j)=0.25*(X^2+T(101,j+1)+2*T(100,j)+T(101,j-1));
             end
             
             for i=2:100
                 T(i,101)=0.25*(X^2+2*T(i,102)+T(i+1,101)+T(i-1,101));
             end
             
             T(101,101)=0.25*(X^2+2*T(101,102)+2*T(102,101));
          %calculation of interior points
             for i=2:100
                 for j=102:125
                     T(i,j)=0.25*(X^2+T(i,j-1)+T(i+1,j)+T(i-1,j)+T(i,j+1));
                 end
             end
             
             for i=77:100
                 for j=126:200
                   T(i,j)=0.25*(X^2+T(i,j-1)+T(i+1,j)+T(i-1,j)+T(i,j+1));
                 end 
             end
        %convective boundary conditions
             for i=77:101
                 T(i,201)=T(i,200)/(1+Bi*X);
             end
             
             for j=126:201
                 T(76,j)=T(77,j)/(1+Bi*X);
             end
             
             for i=2:75
                 T(i,126)=T(i,125)/(1+Bi*X);
             end
             
             for j=101:126
                 T(1,j)=T(2,j)/(1+Bi*X);
             end
            error=max(max(abs(Told-T)));
end
%Implementing symmetry to other quarters
for i=1:101
    j=1:100;
    Told(i,101-j)=Told(i,101+j);
end

for i=1:100
    j=1:201;
    Told(101+i,j)=Told(101-i,j);
end


   