%% CLPT Calculator

clear
clc
close all 

DT=-125; c=0.0125;
Nxyz=[2*1e5; 4.5*1e5; -0.2*1e5];    % Resultant forces given 
Mxyz=[-0.54*1e2; 1e2; 0];           % Resultant torques given 
NM=[Nxyz;Mxyz];

%% Materials of choice and characteristics
% Define the elastic indepenent values for the materials and the structural characteristics of
% each lamina

%moduli(:,:)=[E1       E2       v12      G12        a1          a2       β1      β2] (Pa  C^-1)
moduli(1,:)=[78*1e8  78*1e8    0.32   2954545455   32*1e-6   32*1e-6    0.44    0.44 ];  % CSM isotropic (E1=E2, G=E/2(1+v), a1=a2, β1=β2)
moduli(2,:)=[40*1e9  14*1e9    0.29   4.5*1e9    9.2*1e-6   23.5*1e-6    0       0.44];  %UD ply E-Glass- transeversly isotropic behavior
moduli(3,:)=[40*1e9  14*1e9    0.29   4.5*1e9    9.2*1e-6   23.5*1e-6    0       0.44];  %UD ply E-Glass- transeversly isotropic behavior
moduli(4,:)=[40*1e9  14*1e9    0.29   4.5*1e9    9.2*1e-6   23.5*1e-6    0       0.44];  %UD ply E-Glass- transeversly isotropic behavior
moduli(5,:)=[76*1e9  76*1e9    0.1    5.2*1e9    1.52*1e-6  1.52*1e-6   1e-3   1e-3];  %Woven Roving/Carbon Fibers- orthotropic behavior
moduli(6,:)=[40*1e9  14*1e9    0.29   4.5*1e9    9.2*1e-6   23.5*1e-6    0       0.44];  %UD ply E-Glass- transeversly isotropic behavior
moduli(7,:)=[40*1e9  14*1e9    0.29   4.5*1e9    9.2*1e-6   23.5*1e-6    0       0.44];  %UD ply E-Glass- transeversly isotropic behavior
moduli(8,:)=[40*1e9  14*1e9    0.29   4.5*1e9    9.2*1e-6   23.5*1e-6    0       0.44];  %UD ply E-Glass- transeversly isotropic behavior
moduli(9,:)=[76*1e9  76*1e9    0.1    5.2*1e9    1.52*1e-6  1.52*1e-6   1e-3   1e-3];  %Woven Roving/Carbon Fibers- orthotropic behavior
moduli(10,:)=[78*1e8  78*1e8    0.32   2954545455   32*1e-6   32*1e-6    0.44    0.44 ];  % CSM isotropic (E1=E2, G=E/2(1+v), a1=a2, β1=β2)


%stmoduli(:,:)=[angle thickness] (degrees m)
strmoduli(1,:)=[0   0.6*1e-3];  %CSM
strmoduli(2,:)=[0   0.8*1e-3];  %UD
strmoduli(3,:)=[-45 0.8*1e-3];  %UD
strmoduli(4,:)=[-45 0.8*1e-3];  %UD
strmoduli(5,:)=[ 0  0.4*1e-3];  %WR in plane
strmoduli(6,:)=[ 0  0.8*1e-3];  %UD
strmoduli(7,:)=[-45 0.8*1e-3];  %UD
strmoduli(8,:)=[-45 0.8*1e-3];  %UD
strmoduli(9,:)=[ 0  0.4*1e-3];  %WR in plane
strmoduli(10,:)=[ 0  0.6*1e-3];  %CSM

%strength(:,:)=[X          X'        Y        Y'         Z          Z'       S]
strength(1,:)=[55*1e6   110*1e6   55*1e6   110*1e6     55*1e6   110*1e6    45*1e6];  %CSM
strength(2,:)=[780*1e6  522*1e6   54*1e6   165*1e6     54*1e6   165*1e6    57*1e6];  %UD
strength(3,:)=[780*1e6  522*1e6   54*1e6   165*1e6     54*1e6   165*1e6    57*1e6];  %UD
strength(4,:)=[780*1e6  522*1e6   54*1e6   165*1e6     54*1e6   165*1e6    57*1e6];  %UD
strength(5,:)=[600*1e6  570*1e6   600*1e6  570*1e6    290*1e6   380*1e6    90*1e6];  %WR inplane
strength(6,:)=[780*1e6  522*1e6   54*1e6   165*1e6     54*1e6   165*1e6    57*1e6];  %UD
strength(7,:)=[780*1e6  522*1e6   54*1e6   165*1e6     54*1e6   165*1e6    57*1e6];  %UD
strength(8,:)=[780*1e6  522*1e6   54*1e6   165*1e6     54*1e6   165*1e6    57*1e6];  %UD
strength(9,:)=[600*1e6  570*1e6   600*1e6  570*1e6    290*1e6   380*1e6    90*1e6];  %WR inplane
strength(10,:)=[55*1e6  110*1e6   55*1e6   110*1e6     55*1e6   110*1e6    45*1e6];  %CSM

% ATTENTION!
% When the behaviour is orthotropic (e.g. CSM) X=Y=T=Z and X_=Y_=C=Z_
% (orthotropic->isotropic)
% When the behaviour is transversely isotropic (e.g. UD) Y=Z and Y_=Z_
% and T=S

% Define the number of group repetitions
k=1;

%Calculate the number of the group of laminas
x=length(strmoduli);

%Total number of substrates
n=x*k;

%Define the arrays of angles and sizes
un=zeros(1,n);  %initialization
hn=zeros(1,n);
X=zeros(n,1);
X_=zeros(n,1);
Y=zeros(n,1);
Y_=zeros(n,1);
Z=zeros(n,1);
Z_=zeros(n,1);
S=zeros(n,1);

for i=1:n

    un(1,i)=strmoduli(i,1);
    hn(1,i)=strmoduli(i,2);

    X(i)=strength(i,1);
    X_(i)=strength(i,2);

    Y(i)=strength(i,3);
    Y_(i)=strength(i,4);
    
    Z(i)=strength(i,5);
    Z_(i)=strength(i,6);

    S(i)=strength(i,7);

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


%Effective coefficients
[Ex,Ey,Gxy,vxy,vyx,Exf,Eyf,Gxyt]=laminate_moduli(ABD_t,hn);


%Effective thermal coefficients
[a_eff,b_eff,a_t,b_t,NT,MT,NC,MC,eT]=HTcoeff(un,moduli,n,ABD_t,z,Qnat,hn,DT,c); % a_t=[ax; ay; as] for each lamina and b_t=[bx; by; bs]


% Display the coefficients
disp("                                          ")
disp("Table of coefficients- in plane and flexural:")
table(Ex,Ey,Gxy,vxy,vyx,Exf,Eyf,Gxyt)
disp("                                          ")
disp("Table of hydrothermal coefficients")
table(a_eff,b_eff)


%% Stress distribution for each lamina in the natural coordinate system and the principle coordinate system

Ntotal=Nxyz+NT+NC;  % Total resultant forces
Mtotal=Mxyz+MT+MC;  % Total resultant torques

epso=A_*Ntotal+B_*Mtotal;     % Total real strain of the multilayer structure
kappa=C_*Ntotal+D_*Mtotal;    % Total real curvature of the multilayer structure


[sigma_nat,z_int]=sigma_nat(n,epso,kappa,z,Qnat,a_t,DT,b_t,c);  % Stress distribution of each lamina in the natural coordinate system

sigma_x(:,1)=sigma_nat(1,1:length(sigma_nat));     
sigma_y(:,1)=sigma_nat(2,1:length(sigma_nat));
sigma_s(:,1)=sigma_nat(3,1:length(sigma_nat));

[sigma_prin]= sigma_prin(sigma_nat,un);    % Stress distribution of each lamina in the principle coordinate system

sigma_1(:,1)=sigma_prin(1,1:length(sigma_nat));
sigma_2(:,1)=sigma_prin(2,1:length(sigma_nat));
sigma_6(:,1)=sigma_prin(3,1:length(sigma_nat));

disp("                                          ")
disp("Table of stress distributions in the principle and normal coordinate system accordingly:")
table(z_int,sigma_1,sigma_2,sigma_6,sigma_x,sigma_y,sigma_s)

% Plots of σi-z_int

figure;
hold on;
plot(sigma_x,z_int, 'DisplayName', ' σx ', 'Marker','o', 'LineWidth', 0.8);  
plot(sigma_y,z_int, 'DisplayName', ' σy ', 'Marker','o', 'LineWidth', 0.8);
plot(sigma_s,z_int, 'DisplayName', ' σs ', 'Marker','o', 'LineWidth', 0.8);
grid on;
hold off;
xlabel('σi [Pa]');
ylabel('z [m]');
legend('Location', 'best');

figure;
hold on;
plot(sigma_1,z_int, 'DisplayName', ' σ1 ', 'Marker','o', 'LineWidth', 0.8);  
plot(sigma_2,z_int, 'DisplayName', ' σ2 ', 'Marker','o', 'LineWidth', 0.8);
plot(sigma_6,z_int, 'DisplayName', ' σ6 ', 'Marker','o', 'LineWidth', 0.8);
grid on;
hold off;
xlabel('σi [Pa]');
ylabel('z [m]');
legend('Location', 'best');


%% Failure criteria

% Max stress criterion
[fe]=max_stress(sigma_1,sigma_2,sigma_6,X,X_,Y,Y_,S,n);

% Tsai-Hill criterion
[fTH]=Tsai_Hill(sigma_1,sigma_2,sigma_6,X,X_,Y,Y_,Z,S,n);

% Tsai-Wu criterion
[fTW]=Tsai_Wu(sigma_1,sigma_2,sigma_6,X,X_,Y,Y_,Z,Z_,S,n);

fe_s1(:,1)=fe(1,1:20);
fe_s2(:,1)=fe(2,1:20);
fe_s6(:,1)=fe(3,1:20);

fTH_=transpose(fTH);

fTW_=transpose(fTW);

disp("                                          ")
disp("Table of indices of failure based on Tsai-Hill, EPFS and max stress criteria:")
table(z_int,fTH_,fTW_,fe_s1,fe_s2,fe_s6)


figure;
hold on;
plot(fe_s1,z_int, 'DisplayName', ' σ1max ', 'Marker','*', 'LineWidth', 0.8);  
plot(fe_s2,z_int, 'DisplayName', ' σ2max ', 'Marker','*', 'LineWidth', 0.8);
plot(fe_s6,z_int, 'DisplayName', ' σ6max ', 'Marker','*', 'LineWidth', 0.8);
plot(fTH_,z_int, 'DisplayName', ' Tsai-Hill ', 'Marker','*', 'LineWidth', 0.8);
plot(fTW_,z_int, 'DisplayName', ' EPFS ', 'Marker','*', 'LineWidth', 0.8);
grid on;
hold off;
xlabel('fe');
ylabel('z [m]');
legend('Location', 'best');
