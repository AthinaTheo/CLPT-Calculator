%% CLPT Calculator

clear
clc
close all 

DT=-155; c=0.01125;


%% Materials of choice and characteristics
% Define the elastic indepenent values for the materials and the structural characteristics of
% each lamina

%moduli(:,:)=[E1       E2       v12      G12        a1          a2         β1      β2] (Pa  C^-1)
moduli(1,:)=[181*1e9  10.3*1e9  0.28   7.17*1e9   0.02*1e-6  22.5*1e-6     0     0.44 ];  %  T300/N5208 orthotropic behaviour

%stmoduli(:,:)=[angle thickness] (degrees m)
strmoduli(1,:)=[0   0.6*1e-3];  
strmoduli(2,:)=[45   0.8*1e-3];  
strmoduli(3,:)=[-45 0.8*1e-3];  
strmoduli(4,:)=[90 0.8*1e-3];  

% Define the number of group repetitions
k=1;

%Calculate the number of the group of laminas
X=length(strmoduli);


%Define the arrays of angles and sizes
for i=1:X
    u(1,i)=strmoduli(i,1);
    h(1,i)=strmoduli(i,2);
    
end

%Create the symmetric structures and the total picture of both angles and
%sizes
ut=fliplr(u);
u_values=horzcat(u,ut);                                    
un=repmat(u_values,1,k);  

ht=fliplr(h);
h_values=horzcat(h,ht);
hn=repmat(h_values,1,k);

%Total number of substrates (SYMMETRIC COMPOSITE)
n=2*X*k;

%Define the material properties for each lamina (Made out of the same material here)
for i=1:n
    moduli(i,:)=moduli(1,:);
end

% z array
z=zeros(n+1,1);

ztotal=sum(hn);
z(1)=-ztotal/2;

for j=2:1:n+1
    z(j)=z(j-1)+hn(j-1);
end



% Qprin [stiffness matrix in the principal system]
Qprin=qprin(n,moduli);

% A [extensional stiffness matrix], B [coupling stiffness matrix], D [flexural stiffness matrix]
% Qnat[stiffness matrix in the natural system]
[A,B,D,Qnat]=ABD(n,Qprin,un,hn);

C=transpose(B);
ABD= [A B; C D];


% ABD_t[inversed ABD matrix]
ABD_t=inv(ABD);  

A_=ABD_t(1:3,1:3);
B_=ABD_t(1:3,4:6);
C_=ABD_t(4:6,1:3);
D_=ABD_t(4:6,4:6);



%Effective thermal coefficients
[a_eff,b_eff,a_t,b_t,NT,MT,NC,MC,eT]=HTcoeff(un,moduli,n,ABD_t,z,Qnat,hn,DT,c); % a_t=[ax; ay; as] for each lamina and b_t=[bx; by; bs]


Ntotal=NT+NC;  % Total resultant forces
Mtotal=MT+MC;  % Total resultant torques


epso=A_*Ntotal+B_*Mtotal;     % Total real strain of the multilayer structure
kappa=C_*Ntotal+D_*Mtotal;    % Total real curvature of the multilayer structure


[sigma_nat,z_int]=sigma_nat(n,epso,kappa,z,Qnat,a_t,DT,b_t,c);  % Stress distribution of each lamina in the natural coordinate system

sigma_x(:,1)=sigma_nat(1,1:length(sigma_nat));     
sigma_y(:,1)=sigma_nat(2,1:length(sigma_nat));
sigma_s(:,1)=sigma_nat(3,1:length(sigma_nat));



nexttile
plot(sigma_x,z_int, 'DisplayName', ' σx ');
xlabel('σx');
ylabel('z/ho');
nexttile
plot(sigma_y,z_int,'DisplayName', ' σy ');
xlabel('σy');
ylabel('z/ho');
nexttile
plot(sigma_s,z_int, 'DisplayName', ' σs ' );
xlabel('σs');
ylabel('z/ho');
hold off;



