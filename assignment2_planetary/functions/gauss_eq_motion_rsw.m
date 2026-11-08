function dkepdt = gauss_eq_motion_rsw(t, kep, mu, a_per, e_min)
%GAUSS_EQ_MOTION_RSW_SAFE Gauss planetary equations in the RSW frame (with 1/e regularization).
%
% PROTOTYPE:
%   dkepdt = gauss_eq_motion_rsw(t, kep, mu, a_per, e_min)
%
% PURPOSE:
%   Implements Gauss' planetary equations for the propagation of orbital elements
%   using the components of the perturbing acceleration in the local RSW frame:
%     - R: radial
%     - S: tangential (along-track)
%     - W: normal to the orbital plane (cross-track)
%
% INPUT/OUTPUT STATE:
%   kep = [a; e; i; Omega; omega; M]
%   dkepdt = d/dt[kep]
%
% NUMERICAL NOTES (regularization on e):
%   Some terms in Gauss' equations contain 1/e. To avoid instability when e is very small,
%   e_eff = max(e, e_min) is used ONLY in denominators that include 1/e.
%   (Physical quantities p, h, r, f are computed with the real e.)
%
% INPUT:
%   t     [1]    Time [s]
%   kep   [6x1]  Orbital elements: [a; e; i; Omega; omega; M] [km,-,rad,rad,rad,rad]
%   mu    [1]    Gravitational parameter of the central body [km^3/s^2]
%   a_per [fun]  Handle returning the perturbing acceleration in RSW:
%               a_rsw = a_per(t, kep) = [a_R; a_S; a_W] [km/s^2]
%   e_min [1]    Minimum threshold for regularizing 1/e terms [-]
%
% OUTPUT:
%   dkepdt [6x1] Element derivatives: [da; de; di; dOmega; domega; dM] /dt
%
% DEPENDENCIES:
%   - solve_kepler_E(M,e)  (solution of Kepler's equation for elliptical orbits)
%
% AUTHORS:
%   Simone Frescura
%   Nicola Pifano
%   Luca Puzzi
%   Lorenzo Veronese
%==========================================================================

% ------------------------- State extraction ------------------------------
a     = kep(1);
e     = kep(2);
inc   = kep(3);
Omega = kep(4);
omega = kep(5);
M     = kep(6);

% ------------------ Perturbing acceleration in the RSW frame -------------
a_rsw = a_per(t, kep);
a_r = a_rsw(1);
a_s = a_rsw(2);
a_w = a_rsw(3);

% --------------------- Useful orbital quantities (2BP) -------------------
% p: orbital parameter, h: specific angular momentum
p = a * (1 - e^2);
h = sqrt(mu * p);

% Gauss' equations depend on:
%   - f  (true anomaly)
%   - r  (orbital radius)
% which are obtained from M by solving Kepler: M -> E -> f, then r = p/(1+e cos f)
Mnorm = mod(M + pi, 2*pi) - pi;
E = solve_kepler_E(Mnorm, e);

sinf = (sqrt(1 - e^2) * sin(E)) / (1 - e*cos(E));
cosf = (cos(E) - e) / (1 - e*cos(E));
f = atan2(sinf, cosf);

r = p / (1 + e*cos(f));

% Regularization only for denominators that contain 1/e
e_eff = max(e, e_min);

% -------------------- Gauss planetary equations (RSW) --------------------
da_dt = (2*a^2/h) * ( e*sin(f)*a_r + (p/r)*a_s );

de_dt = (1/h) * ( p*sin(f)*a_r + ( (p + r)*cos(f) + r*e )*a_s );

di_dt = (r*cos(f + omega)/h) * a_w;

dOmega_dt = (r*sin(f + omega)/(h*sin(inc))) * a_w;

domega_dt = (1/(h*e_eff)) * ( -p*cos(f)*a_r + (p + r)*sin(f)*a_s ) ...
         - (r*sin(f + omega)*cos(inc)/(h*sin(inc))) * a_w;

% Mean anomaly derivative (including the perturbing effect)
n = sqrt(mu / a^3);
b = a*sqrt(1 - e^2);

dM_dt = n + (b/(a*h*e_eff)) * ( (p*cos(f) - 2*r*e)*a_r - (p + r)*sin(f)*a_s );

dkepdt = [da_dt; de_dt; di_dt; dOmega_dt; domega_dt; dM_dt];

end

