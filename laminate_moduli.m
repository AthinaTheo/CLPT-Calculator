function [Ex,Ey,Gxy,vxy,vyx,Exf,Eyf,Gxyt]=laminate_moduli(ABD_t,hn)
%% Determine the technical elastic coefficients 

hsum=sum(hn);


Axx_t=ABD_t(1,1); Ayy_t=ABD_t(2,2); Axy_t=ABD_t(1,2); Ass_t=ABD_t(3,3); 

Dxx_t=ABD_t(4,4); Dyy_t=ABD_t(5,5); Dxy_t=ABD_t(4,5); Dss_t=ABD_t(6,6);


%% In plane coefficients
Ex=1/(hsum*Axx_t);  Ey=1/(hsum*Ayy_t);  Gxy=1/(hsum*Ass_t);
vxy=-Axy_t/Axx_t; vyx=-Axy_t/Ayy_t;

%% Flexural coefficients

Exf=12/((hsum^3)*Dxx_t);   Eyf=12/((hsum^3)*Dyy_t);  Gxyt=12/((hsum^3)*Dss_t);

%Row vector of the 5 In-Plane moduli
[moduliplane]=[Ex Ey vxy vyx Gxy];

%Row vector of the 3 Flexural moduli 
[moduliflex]=[Exf Eyf Gxyt];