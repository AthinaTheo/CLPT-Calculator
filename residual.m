function [eR,sigmaR,eM]=residual(A_,B_,C_,D_,NT,MT,z,eT,Qnat,n,epso,kappa)


eot=A_*NT+B_*MT;    %thermal real strain of the middle level of the multilayer in the normal system
kt=C_*NT+D_*MT;     %thermal real curnature of the middle level of the multilayer

for i=1:n
    eps(i,:)=epso+z(i)*kappa;
end

epsT=zeros(n,3);    %thermal real strain of the multilayer eT=[eT1;eT2;eT6]

for i=1:n
    epsT(i,:)=eot+z(i)*kt;
    eR(:,i)=epsT(i,:)-eT(i,:);
    eM(:,i)=eps(i,:)-epsT(i,:);
end

% eT is the stress free thermal strain of each lamina

  % residual strains in each lamina

sigmaR=zeros(n,3);  % residual stresses in each lamina [sigmaxR;sigmayR;sigmasR]

for i=1:n
for j=1:3
    sigmaR(i,:)= Qnat(j,:,i)*eR(:,i);
end
end

      % mechanical real strain of the multilayer

