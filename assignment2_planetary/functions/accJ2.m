function aJ2 = accJ2(r, mu, R, J2)
%ACCJ2 Perturbing acceleration due to the J2 term (oblate gravity field).
%
% PROTOTYPE:
%   aJ2 = accJ2(r, mu, R, J2)
%
% PURPOSE:
%   Computes the perturbing acceleration due to J2 in Cartesian coordinates,
%   assuming a planetocentric equatorial inertial frame (z-axis along the planet
%   rotation axis).
%
% MODEL:
%   Implements the Cartesian components (e.g., eq. 21.16):
%     ax = -(3/2)*mu*J2*(R^2/r^5) * x*(1 - 5*z^2/r^2)
%     ay = -(3/2)*mu*J2*(R^2/r^5) * y*(1 - 5*z^2/r^2)
%     az = -(3/2)*mu*J2*(R^2/r^5) * z*(3 - 5*z^2/r^2)
%
% INPUT:
%   r   [3x1]  Position in the planetocentric equatorial inertial frame [km]
%   mu  [1]    Planet gravitational parameter [km^3/s^2]
%   R   [1]    Mean/equatorial planet radius [km]
%   J2  [1]    Second-degree zonal coefficient [-]
%
% OUTPUT:
%   aJ2 [3x1]  Acceleration due to J2 [km/s^2]
%
% AUTHORS:
%   Simone Frescura
%   Nicola Pifano
%   Luca Puzzi
%   Lorenzo Veronese
%==========================================================================

% Extract the Cartesian components of the position: r = [x; y; z]
x = r(1);
y = r(2);
z = r(3);

% Compute the distance from the planet center r = ||r|| and some recurring quantities:
% - r^2 is needed for the ratio (z^2 / r^2)
% - z^2 appears multiple times in the expressions
rnorm = norm(r);
r2 = rnorm^2;
z2 = z^2;

% Highlight the common factor present in all components of the J2 acceleration:
% K = -(3/2) * mu * J2 * R^2 / r^5
K = -(3/2) * mu * J2 * (R^2) / (rnorm^5);

% Terms in parentheses in the formulas:
% - for ax and ay: (1 - 5*z^2/r^2)
% - for az: (3 - 5*z^2/r^2)
common_xy = (1 - 5*z2/r2);
term_z    = (3 - 5*z2/r2);

% Build the acceleration vector aJ2 = [ax; ay; az]
aJ2 = K * [ x*common_xy;
            y*common_xy;
            z*term_z ];

