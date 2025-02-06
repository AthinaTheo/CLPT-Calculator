%% Function to calculate the Reduced Stiffness Matrix Q of the principal system of each lamina

function [Qprin]= qprin(n,moduli)

for k=1:1:n
    
    E1=moduli(k,1);   %Tensile modulus in the fiber direction 
    E2=moduli(k,2);   %Tensile modulus in the transverse direction 
    v12=moduli(k,3);  %Major Poisson ratio
    G12=moduli(k,4);  %Shear modulus
    v21=v12*E2/E1;    %Minor Poisson ratio  v12>>v21
    
    q=1-(v12*v21);

    Q=[E1/q    (v12*E2)/q   0;
     (v12*E2)/q   E2/q      0;
        0           0       G12];


    for i=1:1:3
        for j=1:1:3
            Qprin(i,j,k)=Q(i,j);
        end
    end
end


