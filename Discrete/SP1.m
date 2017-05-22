%RICHARD VASQUES 

clear all
clc

T=10; % Total length of system
n=100; % # of points
m2=2/((1)^2);
m1=1/1;
%m2=2./((0.5)^2);
c=0.2;
q=1;


h=T/n;
Et=2*m1/m2;
Ea=(1-c)/m1;
D=1/(3*Et);

A=zeros(n+1,n+1);
B=zeros(n+1,1);
SF=zeros(n+1,1);

A(1,1)= 1./2 + 3*D/(2*h);
A(1,2)= -4*D/(2*h);
A(1,3)= D/(2*h);
for i=2:n
    A(i,i-1)= -D/(h^2);
    A(i,i)= Ea+2*D/(h^2);
    A(i,i+1)= -D/(h^2);
    B(i) = q;
end
A(n+1,n-1)= D/(2*h);
A(n+1,n)= -4*D/(2*h);
A(n+1,n+1) = 1./2 + 3*D/(2*h);

SF0=A\B;
SF=1/2*(SF0(2:end)+SF0(1:end-1))
%for i=1:length(SF)
%    SF(i)=1/2*(SF(i)+SF(end+1-i));
%    SF(end+1-i)=SF(i);
%end

p=0+h/2:h:T-h/2;                   % p is line along x-axis.
plot(p,SF,'g'); hold on
%plot(p,SF2,'b');


save SP1.mat
