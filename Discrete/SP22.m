%RICHARD VASQUES 

clear all
clc

T=60; % Total length of system
n=5000; % # of points
m1=1./0.5;
m2=2./(0.5^2);
m3=6./(0.5^3);
m4=24./(0.5^4);
c=1-0.1/900;
q=0.1/900;




l1=3*m4/(10*m2^2)-m3/(3*m1*m2);
b1=m3/(3*m1*m2)-1;

h=T/n;
Et=2*m1/m2;
ea=(1-c)/m1;

pp=ea/Et;
tau=sqrt(2)*sqrt(1+pp+pp^2);
gama=1/6+2*pp/15+tau/8;
beta=tau*gama/(1+4*pp/5);
xi=3*pp*(2*pp-1)/40;
Ea=ea/(1+4*pp/5);
Q=q/(1+4*pp/5);
%Ea=(1-c)*(1-b1*(1-c))/m1;
%Q=q*(1-b1*(1-c));
%D=1/(3*Et);

D=1/(3*Et);


A=zeros(n+1,n+1);
B=zeros(n+1,1);
SF=zeros(n+1,1);
Z=zeros(n+1,1);

A(1,1)= beta + 3*(gama/Et)/(2*h);
A(1,2)= -4*(gama/Et)/(2*h);
A(1,3)= (gama/Et)/(2*h);
B(1)=xi*q/Ea;
for i=2:n
    A(i,i-1)= -D/(h^2);
    A(i,i)= Ea+2*D/(h^2);
    A(i,i+1)= -D/(h^2);
    B(i) = Q;
end
A(n+1,n-1)= (gama/Et)/(2*h);
A(n+1,n)= -4*(gama/Et)/(2*h);
A(n+1,n+1) = beta + 3*(gama/Et)/(2*h);
B(n+1)=xi*q/Ea;

Z=A\B;
SF=(Z+4*q*pp/(5*Ea))./(1+4*pp/5);
SF=1/2*(SF(2:end)+SF(1:end-1));
Z=1/2*(Z(2:end)+Z(1:end-1));
%for i=1:length(SF)
%    SF(i)=1/2*(SF(i)+SF(end+1-i));
%    SF(end+1-i)=SF(i);
%end

p=0+h/2:h:T-h/2;                   % p is line along x-axis.
plot(p,SF,'g'); hold on
%plot(p,SF2,'b');

save SP22.mat
