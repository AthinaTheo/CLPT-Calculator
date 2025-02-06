%% Function to calculate stiffness matrices for CLPT
function [A,B,D,Qnat]=ABD(n,Qprin,un,hn)

Qnat=zeros(3,3,n);
h1=zeros(1,n);
h2=zeros(1,n);
h3=zeros(1,n);
z=zeros(n+1,1);

for i=1:n

    Qnat(1,1,i)=Qprin(1,1,i)*(cosd(un(1,i))^4)+2*(Qprin(1,2,i)+2*Qprin(3,3,i))*(cosd(un(1,i))^2)*(sind(un(1,i))^2)+Qprin(2,2,i)*(sind(un(1,i))^4);
   
    Qnat(1,2,i)=(Qprin(1,1,i)+Qprin(2,2,i)-4*Qprin(3,3,i))*(cosd(un(1,i))^2)*(sind(un(1,i))^2)+Qprin(1,2,i)*((cosd(un(1,i))^4)+(sind(un(1,i))^4));
    
    Qnat(1,3,i)=(Qprin(1,1,i)-Qprin(1,2,i)-2*Qprin(3,3,i))*sind(un(1,i))*(cosd(un(1,i))^3)+(Qprin(1,2,i)-Qprin(2,2,i)+2*Qprin(3,3,i))*(sind(un(1,i))^3)*cosd(un(1,i));
    
    Qnat(2,1,i)=(Qprin(1,1,i)+Qprin(2,2,i)-4*Qprin(3,3,i))*(cosd(un(1,i))^2)*(sind(un(1,i))^2)+Qprin(1,2,i)*((cosd(un(1,i))^4)+(sind(un(1,i))^4));
    
    Qnat(2,2,i)=Qprin(1,1,i)*(sind(un(1,i))^4)+2*(Qprin(1,2,i)+2*Qprin(3,3,i))*(cosd(un(1,i))^2)*(sind(un(1,i))^2)+Qprin(2,2,i)*(cosd(un(1,i))^4);
    
    Qnat(2,3,i)=(Qprin(1,1,i)-Qprin(1,2,i)-2*Qprin(3,3,i))*cosd(un(1,i))*(sind(un(1,i))^3)+(Qprin(1,2,i)-Qprin(2,2,i)+2*Qprin(3,3,i))*(cosd(un(1,i))^3)*sind(un(1,i));
    
    Qnat(3,1,i)=(Qprin(1,1,i)-Qprin(1,2,i)-2*Qprin(3,3,i))*sind(un(1,i))*(cosd(un(1,i))^3)+(Qprin(1,2,i)-Qprin(2,2,i)+2*Qprin(3,3,i))*(sind(un(1,i))^3)*cosd(un(1,i));
    
    Qnat(3,2,i)=(Qprin(1,1,i)-Qprin(1,2,i)-2*Qprin(3,3,i))*cosd(un(1,i))*(sind(un(1,i))^3)+(Qprin(1,2,i)-Qprin(2,2,i)+2*Qprin(3,3,i))*(cosd(un(1,i))^3)*sind(un(1,i));
    
    Qnat(3,3,i)=(Qprin(1,1,i)+Qprin(2,2,i)-2*Qprin(1,2,i)-2*Qprin(3,3,i))*(cosd(un(1,i))^2)*(sind(un(1,i))^2)+Qprin(3,3,i)*((cosd(un(1,i))^4)+(sind(un(1,i))^4));

   
end

ztotal=sum(hn);

z(1)=-ztotal/2;

for j=2:1:n+1
    z(j)=z(j-1)+hn(j-1);
end

A=0;
b=0;
d=0;

for i=1:1:n
    h1(1,i)=(z(i+1)-z(i));
    h2(1,i)=(z(i+1)^2-z(i)^2);
    h3(1,i)=(z(i+1)^3-z(i)^3);

    A=A+(Qnat(:,:,i)*h1(1,i));
    b=b+(Qnat(:,:,i)*h2(1,i));
    d=d+(Qnat(:,:,i)*h3(1,i));
end

B=0.5*b;
D=(1/3)*d;
end