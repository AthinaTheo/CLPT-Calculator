function [fTH]=Tsai_Hill(sigma_1,sigma_2,sigma_6,X,X_,Y,Y_,Z,S,n)

fTH=zeros(1,2*n);
H1=zeros(2*n,1);
H2=zeros(2*n,1);
H12=zeros(2*n,1);
Ath=zeros(2*n,1);

for i=1:n

        if (sigma_1(i*2-1)>0 && sigma_1(1*2)>0)
        H1(i*2-1)=X(i);
        H1(i*2)=X(i);
    else
        H1(i*2-1)=X_(i);
        H1(i*2)=X_(i);
        end

        if (sigma_2(i*2-1)>0 && sigma_2(1*2)>0)
        H2(i*2-1)=Y(i);
        H2(i*2)=Y(i);
    else
        H2(i*2-1)=Y_(i);
        H2(i*2)=Y_(i);
        end

           H12(i*2-1)=(1/(H2(i*2-1)^2))  +(1/(H1(i*2-1)^2)) -(1/(Z(i)^2));

           H12(i*2)=(1/(H2(i*2)^2))  +(1/(H1(i*2)^2)) -(1/(Z(i)^2));
   
 
 Ath(i*2-1)= sigma_1(i*2-1)^2/(H1(i*2-1)^2) - H12(i*2-1)*sigma_1(i*2-1)*sigma_2(i*2-1)  ...
             + sigma_2(i*2-1)^2/(H2(i*2-1)^2)  + (sigma_6(i*2-1)^2/S(i)^2);

 Ath(i*2)= sigma_1(i*2)^2/(H1(i*2)^2) - H12(i*2)*sigma_1(i*2)*sigma_2(i*2)  ...
             + sigma_2(i*2)^2/(H2(i*2)^2)  + (sigma_6(i*2)^2/S(i)^2);


 fTH(1,i*2-1)=sqrt(Ath(i*2-1));
 fTH(1,i*2)=sqrt(Ath(i*2));
end
