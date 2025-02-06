function [fe]=max_stress(sigma_1,sigma_2,sigma_6,X,X_,Y,Y_,S,n)

fe=zeros(3,2*n);

for i=1:n

        if (sigma_1(i*2-1)>0 && sigma_1(1*2)>0)
        fe(1,i*2-1)=sigma_1(i*2-1)/X(i);
        fe(1,i*2)=sigma_1(i*2)/X(i);
    else
        fe(1,i*2-1)=abs(sigma_1(i*2-1))/X_(i);
        fe(1,i*2)=abs(sigma_1(i*2))/X_(i);
        end

        if (sigma_2(i*2-1)>0 && sigma_2(1*2)>0)
        fe(2,i*2-1)=sigma_2(i*2-1)/Y(i);
        fe(2,i*2)=sigma_2(i*2)/Y(i);
    else
        fe(2,i*2-1)=abs(sigma_2(i*2-1))/Y_(i);
        fe(2,i*2)=abs(sigma_2(i*2))/Y_(i);
        end

        if (sigma_6(i*2-1)>0 && sigma_6(1*2)>0)
        fe(3,i*2-1)=sigma_6(i*2-1)/S(i);
        fe(3,i*2)=sigma_6(i*2)/S(i);
    else
        fe(3,i*2-1)=abs(sigma_6(i*2-1))/S(i);
        fe(3,i*2)=abs(sigma_6(i*2))/S(i);
        end        
end
