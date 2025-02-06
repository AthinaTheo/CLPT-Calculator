function [fTW]=Tsai_Wu(sigma_1,sigma_2,sigma_6,X,X_,Y,Y_,Z,Z_,S,n)

fTW=zeros(1,2*n);
b=zeros(2*n,1);
c=zeros(2*n,1);
F11=zeros(n,1);
F22=zeros(n,1);
F33=zeros(n,1);
F66=zeros(n,1);
F12=zeros(n,1);
F1=zeros(n,1);
F2=zeros(n,1);

for i=1:n

    F11(i)=1/(X(i)*X_(i));
    F22(i)=1/(Y(i)*Y_(i));
    F33(i)=1/(Z(i)*Z_(i));
    F66(i)=1/S(i)^2;
    F12(i)=0.5*(F33(i)-F11(i)-F22(i));
    F1(i)=(1/X(i))-(1/X_(i));
    F2(i)=(1/Y(i))-(1/Y_(i));


    b(2*i-1)=F1(i)*sigma_1(2*i-1)+F2(i)*sigma_2(2*i-1);
    b(2*i)=F1(i)*sigma_1(2*i)+F2(i)*sigma_2(2*i);

    c(2*i-1)=F11(i)*sigma_1(2*i-1)^2+F22(i)*sigma_2(2*i-1)^2+F66(i)*sigma_6(2*i-1)^2+2*F12(i)*sigma_1(2*i-1)*sigma_2(2*i-1);
    c(2*i)=F11(i)*sigma_1(2*i)^2+F22(i)*sigma_2(2*i)^2+F66(i)*sigma_6(2*i)^2+2*F12(i)*sigma_1(2*i)*sigma_2(2*i);

    fTW(1,2*i-1)=0.5*(b(2*i-1)+sqrt(b(2*i-1)^2+4*c(2*i-1)));
    fTW(1,2*i)=0.5*(b(2*i)+sqrt(b(2*i)^2+4*c(2*i)));
end
    
