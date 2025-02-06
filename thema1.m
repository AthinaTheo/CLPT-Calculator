%% CLPT Calculator

clear
clc
close all 

%% Materials of choice and characteristics
% Define the elastic indepenent values for the materials and the structural characteristics of
% each lamina

%moduli(:,:)=[E1       E2       v12      G12 ] (Pa  C^-1)
moduli(1,:)=[78*1e8  78*1e8    0.32   2954545455];% CSM isotropic (E1=E2, G=E/2(1+v))
moduli(2,:)=[78*1e8  78*1e8    0.32   2954545455];  % CSM isotropic
moduli(3,:)=[38*1e9  15*1e9    0.25   5.5*1e9];  %UD ply E-Glass- transeversly isotropic behavior
moduli(4,:)=[38*1e9  15*1e9    0.25   5.5*1e9];  %UD ply E-Glass- transeversly isotropic behavior
moduli(5,:)=[38*1e9  15*1e9    0.25   5.5*1e9];  %UD ply E-Glass- transeversly isotropic behavior
moduli(6,:)=[96*1e9  96*1e9    0.03    7.2*1e9];  %Woven Roving/Carbon Fibers- orthotropic behavior

%stmoduli(:,:)=[angle thickness] (degrees m)
strmoduli(1,:)=[0   0.6*1e-3];  %CSM
strmoduli(2,:)=[0   0.6*1e-3];  %CSM
strmoduli(3,:)=[0   0.8*1e-3];  %UD
strmoduli(4,:)=[-45 0.8*1e-3];  %UD
strmoduli(5,:)=[-45 0.8*1e-3];  %UD
strmoduli(6,:)=[ 0  0.4*1e-3];  %WR in plane

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

Qprin=qprin(n,moduli);

[A,B,D,Qnat]=ABD(n,Qprin,un,hn);

[Ex,Ey,Gxy,vxy,vyx,Exf,Eyf,Gxyt]=laminate_moduli(A,D,hn);

%% Determination of size of the structure-plate
a=2.7; b=3.9; R=a/b;

Dxx=D(1,1); Dxy=D(1,2); Dxs=D(1,3); 
Dyy=D(2,2); Dys=D(2,3); Dss=D(3,3); 

pmn1=1;
m1=2; 
n1=1; 

Dmn1=Dxx*m1^4+2*(Dxy+2*Dss)*(m1*n1*R)^2+Dyy*(n1*R)^4;


wmn1=(a^4/pi^4)*pmn1/Dmn1;


x=0:0.027:2.67;
y=0:0.039:3.9;

%Initialize the matrices

Mx=zeros(length(x),length(y)); My=zeros(length(x),length(y)); Ms=zeros(length(x),length(y)); 


Qx=zeros(length(x),length(y)); Qy=zeros(length(x),length(y));


kx=zeros(length(x),length(y)); ky=zeros(length(x),length(y)); ks=zeros(length(x),length(y));


d2w_dx2=zeros(length(x),length(y)); d2w_dy2=zeros(length(x),length(y)); d2w_dxdy=zeros(length(x),length(y));


d3w_dx3=zeros(length(x),length(y)); d3w_dy3=zeros(length(x),length(y)); d3w_dy2dx=zeros(length(x),length(y)); d3w_dx2dy=zeros(length(x),length(y));


dMx_dx=zeros(length(x),length(y)); dMx_dy=zeros(length(x),length(y)); dMy_dx=zeros(length(x),length(y)); 
dMy_dy=zeros(length(x),length(y)); dMs_dx=zeros(length(x),length(y)); dMs_dy=zeros(length(x),length(y));
                               
w=zeros(length(x),length(y));

p=zeros(length(x),length(y));


%Boundary conditions
Mx(1,1)=0; Mx(end,1)=0; My(1,1)=0; My(end,1)=0;
po=1;

%Auxiliary variables
Co=a^4/pi^4*(pmn1/Dmn1); C1=m1*pi/a; L1=n1*pi/b; 


%Loop to calculate the resultant forces and torques 
for i=1:length(x)

  for j=1:length(y)
%Pressure p
p(i,j)=po*(pmn1*sin(C1*x(i))*sin(L1*y(j)));

%Function that determines how bending affects each area
w(i,j)=po*Co*sin(C1*x(i))*sin(L1*y(j));

  end
  %Same procedure for i index
p(i,j)=po*(pmn1*sin(C1*x(i))*sin(L1*y(j)));

w(i,j)=Co*sin(C1*x(i))*sin(L1*y(j));


end



%% Distribution of stresses 

%Definition of maximum plundge and its position
wmin=min(w,[],"all");
[minRow,minCol]=find(w==wmin);
