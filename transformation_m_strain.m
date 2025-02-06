function [Te]=transformation_m_strain(un)

c=cosd(un);  % angle in degrees
s=sind(un);  % angle in degrees

Te=[c^2 s^2 c*s;
    s^2 c^2 -c*s;
    -2*c*s 2*c*s c^2-s^2];
