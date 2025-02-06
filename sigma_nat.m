function [sigma_nat,z_int]=sigma_nat(n,epso,kappa,z,Qnat,a_t,DT,b_t,c)


z_int=repelem(z,2);               % Array of length 2*n with the edge location of each ply
z_int=z_int(2:end-1);
sigma_nat=zeros(3,length(z_int)); % Matrix of size 3x2n with the stress state (natural coordinate system) in the plies at the top and bottom of each ply. 

eps=zeros(2*n,3);                 % Total real strain of the multilayer of each lamina + thermal strain of each lamina eiT=a_t*DT

for i=1:n
    eps=[epso+z_int(i*2-1)*kappa-a_t(:,i)*DT-b_t(:,i)*c, epso+z_int(i*2)*kappa-a_t(:,i)*DT-b_t(:,i)*c];
    sigma_nat(:,i*2-1)= Qnat(:,:,i)*eps(:,1);
    sigma_nat(:,i*2)=Qnat(:,:,i)*eps(:,2);
end





