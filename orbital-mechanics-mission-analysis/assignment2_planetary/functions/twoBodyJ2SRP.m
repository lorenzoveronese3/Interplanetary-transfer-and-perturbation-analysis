function dy = twoBodyJ2SRP(t, y, mu, R, J2, mjd2000_0, ibody, Cr, AoverM)
%TWOBODYJ2SRP ODE system in Cartesian coordinates: 2BP + J2 + SRP perturbations.
%
% PROTOTYPE:
%   dy = twoBodyJ2SRP(t, y, mu, R, J2, mjd2000_0, ibody, Cr, AoverM)
%
% PURPOSE:
%   Returns the derivative of the Cartesian state y = [r; v] including:
%   central gravity (2BP), J2 perturbation, and solar radiation pressure (SRP).
%
% DYNAMICAL MODEL:
%   r_dot = v
%   v_dot = a_2BP(r) + a_J2(r) + a_SRP(t,r)
%
% INPUT:
%   t         [1]    Time from the initial epoch [s]
%   y         [6x1]  Cartesian state [r; v] with r [km], v [km/s]
%   mu        [1]    Gravitational parameter [km^3/s^2]
%   R         [1]    Mean planet radius [km]
%   J2        [1]    J2 coefficient [-]
%   mjd2000_0 [1]    Initial epoch in MJD2000 [days]
%   ibody     [1]    Planet ID [-]
%   Cr        [1]    SRP reflectivity coefficient [-]
%   AoverM    [1]    SRP area-to-mass ratio [m^2/kg]
%
% OUTPUT:
%   dy        [6x1]  State derivative [rdot; vdot] with vdot in [km/s^2]
%
% AUTHORS:
%   Simone Frescura
%   Nicola Pifano
%   Luca Puzzi
%   Lorenzo Veronese
%==========================================================================

% State: position and velocity
r = y(1:3);     % [km]
v = y(4:6);     % [km/s]

% Accelerations: 2BP (central) + J2 (oblateness) + SRP (cannonball)
a2BP = -mu * r / (norm(r)^3);
aJ2  = accJ2(r, mu, R, J2);

% SRP: Sun direction obtained from uplanet; frame conversion handled in accSRP
aSRP = accSRP(t, r, mjd2000_0, ibody, Cr, AoverM);

% ODE system: rdot = v, vdot = sum of accelerations
dy = [v; (a2BP + aJ2 + aSRP)];

end



