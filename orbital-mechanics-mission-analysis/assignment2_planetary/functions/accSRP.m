function aSRP = accSRP(t, r, mjd2000_0, ibody, Cr, AoverM)
%ACCSRP Acceleration due to Solar Radiation Pressure (SRP) - cannonball model.
%
% PROTOTYPE:
%   aSRP = accSRP(t, r, mjd2000_0, ibody, Cr, AoverM)
%
% PURPOSE:
%   Computes the SRP acceleration in the planetocentric equatorial inertial frame.
%   The Sun direction is obtained from uplanet (ecliptic frame) and is
%   rotated into the equatorial frame via the planet axial tilt.
%
% INPUT:
%   t         [1]    Time from the initial epoch [s]
%   r         [3x1]  Satellite position (planetocentric inertial) [km]
%   mjd2000_0 [1]    Initial epoch in MJD2000 [days]
%   ibody     [1]    Planet ID (consistent with uplanet) [-]
%   Cr        [1]    Reflectivity coefficient [-]
%   AoverM    [1]    Area-to-mass ratio [m^2/kg]
%
% OUTPUT:
%   aSRP      [3x1]  SRP acceleration [km/s^2]
%
% AUTHORS:
%   Simone Frescura
%   Nicola Pifano
%   Luca Puzzi
%   Lorenzo Veronese
%==========================================================================

% Absolute time in MJD2000
mjd2000 = mjd2000_0 + t/86400;

% Heliocentric state of the planet (ecliptic frame)
[kepP, muSun] = uplanet(mjd2000, ibody);

% rP_S = Sun -> Planet (ecliptic) [km]
[rP_S_ecl, ~] = kep2car(kepP, muSun);

% rS_P = Planet -> Sun (ecliptic) [km]
rS_P_ecl = -rP_S_ecl;

% Ecliptic -> planetary equatorial rotation via axial tilt
eps = deg2rad(astroConstants(60 + ibody));
Rx = [ 1      0        0;
       0   cos(eps) -sin(eps);
       0   sin(eps)  cos(eps) ];
rS_P = Rx * rS_P_ecl;

% Unit vector satellite -> Sun:
% rS_P is (Planet -> Sun) and r is (Planet -> Satellite), so (rS_P - r) = (Satellite -> Sun)
d   = rS_P - r;           % [km]
e_s = d / norm(d);

% Irradiance at the planet [W/m^2] and speed of light [m/s]
W     = astroConstants(70 + ibody);
c_mps = astroConstants(5) * 1000;

% Cannonball SRP (push away from the Sun)
% e_s points toward the Sun; the "-" sign makes the acceleration directed away from the Sun (photon pressure)
a_mps = -Cr * (W / c_mps) * AoverM * e_s;  % [m/s^2]
aSRP  = a_mps / 1000;                      % [km/s^2]

end







