function a_rsw = accPert_J2_SRP_RSW(t, kep, pars)
%ACCPERT_J2_SRP_RSW Computes the perturbing acceleration (J2 + SRP) in the RSW frame.
%
% PROTOTYPE:
%   a_rsw = accPert_J2_SRP_RSW(t, kep, pars)
%
% PURPOSE:
%   Computes the total perturbing acceleration aPert = aJ2 + aSRP and expresses it
%   in the local RSW reference frame associated with the orbit:
%     - R: radial direction (from the planet to the satellite)
%     - S: tangential direction (along-track)
%     - W: normal direction to the orbital plane (cross-track)
%
% DESCRIPTION (what it does and why):
%   The input kep is given with mean anomaly M, which does not directly provide the position.
%   To obtain r and v, the following steps are needed:
%     1) M -> E  (Kepler's equation)
%     2) E -> θ  (true anomaly, to build the "classic" Keplerian state)
%     3) [a,e,i,Ω,ω,θ] -> r,v  (Keplerian->Cartesian conversion)
%   Once r,v are known:
%     4) compute aJ2(r) and aSRP(t,r) in Cartesian coordinates
%     5) build the RSW triad from r and v and project aPert onto [R,S,W]
%
% INPUT:
%   t    [1]    Time from the initial epoch mjd2000_0 [s]
%   kep  [6x1]  Elements: [a; e; i; Omega; omega; M] [km,-,rad,rad,rad,rad]
%   pars [1]    Parameter structure with fields:
%              mu, R, J2, mjd2000_0, ibody, Cr, AoverM
%
% OUTPUT:
%   a_rsw [3x1] Perturbing acceleration in RSW [km/s^2] = [a_R; a_S; a_W]
%
% AUTHORS:
%   Simone Frescura
%   Nicola Pifano
%   Luca Puzzi
%   Lorenzo Veronese
%==========================================================================


% 1) Orbital elements (with M as mean anomaly)
a     = kep(1);
e     = kep(2);
inc   = kep(3);
Omega = kep(4);
omega = kep(5);
M     = kep(6);

% 2) M -> E by solving Kepler's equation: M = E - e*sin(E)
Mnorm = mod(M + pi, 2*pi) - pi;
E = solve_kepler_E(Mnorm, e);

% 3) E -> theta (true anomaly) with a robust formula (atan2)
sinf  = (sqrt(1 - e^2) * sin(E)) / (1 - e*cos(E));
cosf  = (cos(E) - e) / (1 - e*cos(E));
theta = atan2(sinf, cosf);
if theta < 0, theta = theta + 2*pi; end

% 4) Reconstruct r,v from the classic Keplerian state [a,e,i,Ω,ω,θ]
kep_theta = [a; e; inc; Omega; omega; theta];
[r, v] = kep2car(kep_theta, pars.mu);

% 5) Total perturbing acceleration in Cartesian coordinates
aJ2  = accJ2(r, pars.mu, pars.R, pars.J2);
aSRP = accSRP(t, r, pars.mjd2000_0, pars.ibody, pars.Cr, pars.AoverM);
aP = aJ2 + aSRP;

% 6) Build the RSW triad from r and v:
%    R = r/||r||, W = (r x v)/||r x v||, S = W x R
rhat = r / norm(r);
hvec = cross(r, v);
what = hvec / norm(hvec);
shat = cross(what, rhat);

% 7) Project aP onto the R, S, W directions
a_r = dot(aP, rhat);
a_s = dot(aP, shat);
a_w = dot(aP, what);

a_rsw = [a_r; a_s; a_w];
end









