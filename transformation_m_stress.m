function [Ts]=transformation_m_stress(un)

c=cosd(un);  % angle in degrees
s=sind(un);  % angle in degrees

Ts=[c^2 s^2 2*c*s;
    s^2 c^2 -2*c*s;
    -c*s c*s c^2-s^2];
