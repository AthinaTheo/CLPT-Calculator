function [a_eff,b_eff,a_t,b_t,NT,MT,NC,MC,eT]=HTcoeff(un,moduli,n,ABD_t,z,Qnat,hn,DT,c)

for i=1:1:n
a(1,i)=moduli(i,5);
a(2,i)=moduli(i,6);
b(1,i)=moduli(i,7);
b(2,i)=moduli(i,8);

ax(1,i)=a(1,i)*cosd(un(1,i))^2+a(2,i)*sind(un(1,i))^2;
ay(1,i)=a(1,i)*sind(un(1,i))^2+a(2,i)*cosd(un(1,i))^2;
as(1,i)=2*sind(un(1,i))*cosd(un(1,i))*(a(1,i)-a(2,i));

bx(1,i)=b(1,i)*cosd(un(1,i))^2+b(2,i)*sind(un(1,i))^2;
by(1,i)=b(1,i)*sind(un(1,i))^2+b(2,i)*cosd(un(1,i))^2;
bs(1,i)=2*sind(un(1,i))*cosd(un(1,i))*(b(1,i)-b(2,i));
end

a_t=[ax; ay; as];
b_t=[bx; by; bs];

Qak1=0;
Qak2=0;
Qbk1=0;
Qbk2=0;
NT=0;
MT=0;
NC=0;
MC=0;
eT=zeros(n,3);

for i=1:n
    Qak1=Qak1+(Qnat(:,:,i)*a_t(:,i)*(z(i+1)-z(i)));
    Qak2=Qak2+(Qnat(:,:,i)*a_t(:,i)*(z(i+1)^2/2-z(i)^2/2));
    Qbk1=Qbk1+(Qnat(:,:,i)*b_t(:,i)*(z(i+1)-z(i)));
    Qbk2=Qbk2+(Qnat(:,:,i)*b_t(:,i)*(z(i+1)^2/2-z(i)^2/2));

    NT=NT+(Qnat(:,:,i)*a_t(:,i)*(z(i+1)-z(i)))*DT;
    MT=MT+(Qnat(:,:,i)*a_t(:,i)*(z(i+1)^2/2-z(i)^2/2))*DT;
    NC=NC+(Qnat(:,:,i)*b_t(:,i)*(z(i+1)-z(i)))*c;
    MC=MC+(Qnat(:,:,i)*b_t(:,i)*(z(i+1)^2/2-z(i)^2/2))*c;


    eT(i,:)=a_t(:,i)*DT;

end


aox=ABD_t(1,1)*Qak1(1)+ABD_t(1,2)*Qak1(2)+ABD_t(1,3)*Qak1(3)  +ABD_t(4,1)*Qak2(1)+ABD_t(5,1)*Qak2(2)+ABD_t(6,1)*Qak2(3);
aoy=ABD_t(2,1)*Qak1(1)+ABD_t(2,2)*Qak1(2)+ABD_t(2,3)*Qak1(3)  +ABD_t(4,2)*Qak2(1)+ABD_t(5,2)*Qak2(2)+ABD_t(6,2)*Qak2(3);
aos=ABD_t(3,1)*Qak1(1)+ABD_t(3,2)*Qak1(2)+ABD_t(3,3)*Qak1(3)  +ABD_t(4,3)*Qak2(1)+ABD_t(5,3)*Qak2(2)+ABD_t(6,3)*Qak2(3);

a_eff=[aox;aoy;aos];

box=ABD_t(1,1)*Qbk1(1)+ABD_t(1,2)*Qbk1(2)+ABD_t(1,3)*Qbk1(3)  +ABD_t(4,1)*Qbk2(1)+ABD_t(5,1)*Qbk2(2)+ABD_t(6,1)*Qbk2(3);
boy=ABD_t(2,1)*Qbk1(1)+ABD_t(2,2)*Qbk1(2)+ABD_t(2,3)*Qbk1(3)  +ABD_t(4,2)*Qbk2(1)+ABD_t(5,2)*Qbk2(2)+ABD_t(6,2)*Qbk2(3);
bos=ABD_t(3,1)*Qbk1(1)+ABD_t(3,2)*Qbk1(2)+ABD_t(3,3)*Qbk1(3)  +ABD_t(4,3)*Qbk2(1)+ABD_t(5,3)*Qbk2(2)+ABD_t(6,3)*Qbk2(3);

b_eff=[box;boy;bos];

