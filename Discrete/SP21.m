%RICHARD VASQUES 

clear all
clc

T=60; % Total length of system
n=120; % # of points
m1=1/0.5;
m2=2/(0.5)^2;
m3=6/(0.5)^3;
m4=24/(0.5)^4;
c=1-0.4995/0.5;
q = 1;


l1=3*m4/(10*m2^2)-m3/(3*m1*m2);
b1=m3/(3*m1*m2)-1;

h=T/n;
Et=2*m1/m2;
Ea=((1-c)/m1)*(1-b1*(1-c))/(1+l1*(1-c));
Q=q*(1-b1*(1-c))/(1+l1*(1-c));
D=1/(3*Et);

A=zeros(n+1,n+1);
B=zeros(n+1,1);
SF=zeros(n+1,1);
Z=zeros(n+1,1);


A(1,1)= 1./2 + 3*D/(2*h);
A(1,2)= -4*D/(2*h);
A(1,3)= D/(2*h);
for i=2:n
    A(i,i-1)= -D/(h^2);
    A(i,i)= Ea+2*D/(h^2);
    A(i,i+1)= -D/(h^2);
    if i*h>(T/2-0.5)
        if i*h<(T/2+0.5)
            B(i) = Q;
        end
    end
end
A(n+1,n-1)= D/(2*h);
A(n+1,n)= -4*D/(2*h);
A(n+1,n+1) = 1./2 + 3*D/(2*h);

Z=A\B;

for i=2:n
    if i*h>(T/2-0.5)
        if i*h<(T/2+0.5)
            Z(i)=(Z(i)+l1*m1*q);
        end
    end
end
SF(:)=(Z(:))./(1+l1*(1-c));
SF=1/2*(SF(2:end)+SF(1:end-1))
%for i=1:length(SF)
%    SF(i)=1/2*(SF(i)+SF(end+1-i));
%    SF(end+1-i)=SF(i);
%end

p=0+h/2:h:T-h/2;                   % p is line along x-axis.
plot(p,SF,'b'); hold on
%plot(p,SF2,'k');

save SP2b.mat