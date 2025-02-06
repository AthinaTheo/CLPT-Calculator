function [sigma_prin]= sigma_prin(sigma_nat,un)

sigma_prin= zeros(size(sigma_nat));

for i=1:length(un)
    T=transformation_m_stress(un(1,i));
    sigma_prin(:,i*2-1)=T*sigma_nat(:,i*2-1);
    sigma_prin(:,i*2)=T*sigma_nat(:,i*2);
end