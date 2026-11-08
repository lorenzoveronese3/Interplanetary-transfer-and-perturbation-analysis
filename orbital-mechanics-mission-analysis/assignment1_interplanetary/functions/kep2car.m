function [r_vec,v_vec] = kep2car(a,e,i,OMEGA,omega,theta,mu)
% This function convert keplerian elements into cartesian coordinates.
% INPUT: 
%       - a[1]: semi-major axis [km]
%       - e[1]: eccentricity (scalar) [-]
%       - i[1]: inclination [rad]
%       - OMEGA[1]: right ascension of the ascending node [rad]
%       - omega[1]: argument of the pericentre [rad]
%       - theta[1]: true anomaly [rad]
%       - mu[1]: gravitational constant [km^3/s^2]
%
% OUTPUT:
%        - r_vec [3]: position vector [km]
%        - v_vec [3]: velocity vector [km/s]
%
% CONTRIBUTORS:
% -Simone Frescura
% -Lorenzo Veronese
% -Nicola Pifano
% -Luca Puzzi
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

% Definition of r end v in perifocal frame
p=a*(1-e^2);
r=p/(1+e*cos(theta));
r_per=r*[cos(theta); sin(theta); 0];
v_per= sqrt(mu/p) * [-sin(theta); e+cos(theta); 0];

% Rotation Matrices

R_3om = [cos(omega) sin(omega) 0;   
        -sin(omega) cos(omega) 0;
         0           0         1];
R_1i = [ 1    0      0;
         0 cos(i) sin(i);
         0 -sin(i) cos(i)];
R_3OM = [cos(OMEGA) sin(OMEGA) 0;
        -sin(OMEGA) cos(OMEGA) 0;
         0           0         1];

% Rotation composition from perifocal to inertial
Q=(R_3om*R_1i*R_3OM)';

% output
r_vec = Q*r_per;
v_vec = Q*v_per;

end