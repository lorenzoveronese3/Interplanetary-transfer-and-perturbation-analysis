function [a,e,i,OMEGA,omega,theta] = car2kep(r_vec,v_vec,mu)
% This function compute the keplerian elements from the state vector in ECIEq
% Input:
%        - r_vec[3]: position vector [km]
%        - v_vec[3]: velocity vector [km/s]
%        - mu[1]: gravitational constant of the planet [km^3/s^2]
%
% Output:
%        - a[1]: semi-major axis [km]
%        - e[1]: eccentricity (scalar) [-]
%        - i[1]: inclination [rad]
%        - OMEGA[1]: right ascension of the ascending node [rad]
%        - omega[1]: argument of the pericentre [rad]
%        - theta[1]: true anomaly [rad]
%
% CONTRIBUTORS:
% -Simone Frescura
% -Lorenzo Veronese
% -Nicola Pifano
% -Luca Puzzi
%
% -------------------------------------------------------------------------

% Angular momentum h
h_vec=cross(r_vec,v_vec);
h=norm(h_vec);

% Inclination i
i=acos(h_vec(end)/h);  %[rad]

% Line of nodes N and Right Ascension of Ascending Node OMEGA
if i==0
    N_vers=[1 0 0]';
    OMEGA=0;
else
    N_vec=cross([0 0 1]',h_vec);
    N=norm(N_vec);
    N_vers=N_vec/N;
    if N_vers(2)>=0
        OMEGA=acos(N_vers(1));
    else
        OMEGA=2*pi-acos(N_vers(1));
    end
end

% Eccentricity e
e_vec=cross(v_vec,h_vec)/mu - r_vec/norm(r_vec);
e=norm(e_vec);

% Argument of the pericentre omega
if e<1e-5
    omega=0;
else
    if e_vec(3)>=0
        omega=acos(dot(N_vers,e_vec/e));
    else
        omega=2*pi - acos(dot(N_vers,e_vec/e));
    end
end

% True anomaly theta
v_rad=dot(v_vec,r_vec/norm(r_vec));  % radial velocity
if v_rad>=0
    theta=acos(dot(e_vec,r_vec)/(e*norm(r_vec)));
else
    theta=2*pi-acos(dot(e_vec,r_vec)/(e*norm(r_vec)));
end

% Semi-major axis a
a=mu*norm(r_vec)/(2*mu-norm(r_vec)*dot(v_vec,v_vec));

end