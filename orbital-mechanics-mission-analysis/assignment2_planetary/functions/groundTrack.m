function [alpha, delta, lon, lat] = groundTrack(state0, theta0, tVec, mu, omegaE, varargin)
%GROUNDTRACK Ground track computation in the 2BP model.
%
% PROTOTYPE:
%   [alpha, delta, lon, lat] = groundTrack(state0, theta0, tVec, mu, omegaE, ...)
%
% PURPOSE:
%   1) Propagates the orbit with the two-body model (2BP) in Cartesian coordinates.
%   2) Converts the inertial position r(t) into right ascension alpha(t) and declination delta(t).
%   3) Converts (alpha,delta) into (lon,lat) in the planet-fixed frame using:
%        lon = alpha - theta(t),  lat = delta
%      where theta(t) is the planet prime meridian angle.
%
% DYNAMICAL MODEL (2BP):
%   r_dot = v
%   v_dot = -mu * r / ||r||^3
%
% INPUT:
%   state0 [6x1]  Initial state:
%                - 'cart': [r0; v0] with r0 [km], v0 [km/s]
%                - 'kep' : [a; e; i; RAAN; omega; theta] (angles in [rad])
%   theta0 [1]    Planet prime meridian angle at the initial time [rad]
%   tVec   [Nx1]  Output times [s]
%   mu     [1]    Gravitational parameter of the central body [km^3/s^2]
%   omegaE [1]    Angular velocity used for theta(t) in the ground track [rad/s]
%                (follows the course convention; can be < 0 for retrograde rotation)
%
%   Name-Value options:
%     'InputType'  : 'cart' (default) | 'kep'
%     'OutputUnit' : 'deg'  (default) | 'rad'
%     'RelTol'     : default 1e-11
%     'AbsTol'     : default 1e-12
%     't0'         : reference for theta(t), default tVec(1)
%
% OUTPUT:
%   alpha, delta [Nx1] Right ascension and declination (inertial frame) [deg or rad]
%   lon, lat     [Nx1] Longitude and latitude (planet-fixed frame)      [deg or rad]
%
% DEPENDENCIES:
%   - twoBody2BP(t,y,mu)
%   - kep2car(kep,mu)
%
% AUTHORS:
%   Simone Frescura
%   Nicola Pifano
%   Luca Puzzi
%   Lorenzo Veronese
%==========================================================================

% --------------------------- Name-Value options --------------------------
opt.InputType  = 'cart';
opt.OutputUnit = 'deg';
opt.RelTol     = 1e-11;
opt.AbsTol     = 1e-12;
opt.t0         = tVec(1);

for k = 1:2:numel(varargin)
    opt.(varargin{k}) = varargin{k+1};
end

tVec = tVec(:);

% --------------------- Initial state in Cartesian ------------------------
% If requested, convert from Keplerian to Cartesian for the 2BP propagation.
if strcmpi(opt.InputType,'kep')
    [r0, v0] = kep2car(state0(:), mu);
    y0 = [r0; v0];
else
    y0 = state0(:);
end

% ------------------------- 1) 2BP propagation ----------------------------
odefun = @(t,y) twoBody2BP(t, y, mu);
odeOpt = odeset('RelTol', opt.RelTol, 'AbsTol', opt.AbsTol);
[~, Y] = ode113(odefun, tVec, y0, odeOpt);

% ----------------- 2) r(t) -> inertial (alpha,delta) ---------------------
% (alpha, delta) are the spherical coordinates of the position in the inertial frame:
%   alpha = atan2(y,x)  (right ascension)
%   delta = asin(z/||r||) (declination)
x = Y(:,1); y = Y(:,2); z = Y(:,3);
rnorm = sqrt(x.^2 + y.^2 + z.^2);

delta_rad = asin(z./rnorm);
alpha_rad = mod(atan2(y,x), 2*pi);

% -------------- 3) inertial -> planet-fixed: (lon,lat) -------------------
% theta(t) is the prime meridian angle of the planet over time:
%   theta(t) = theta0 + omegaE*(t - t0)
% Course ground track convention:
%   lon = alpha - theta(t),  lat = delta
dt = tVec - opt.t0;
theta = mod(theta0 + omegaE.*dt, 2*pi);

lon_rad = alpha_rad - theta;
lon_rad = mod(lon_rad + pi, 2*pi) - pi;   % wrap in [-pi, +pi] to avoid jumps in plots
lat_rad = delta_rad;

% --------------------------- Output units --------------------------------
if strcmpi(opt.OutputUnit,'rad')
    alpha = alpha_rad;
    delta = delta_rad;
    lon   = lon_rad;
    lat   = lat_rad;
else
    alpha = rad2deg(alpha_rad);
    delta = rad2deg(delta_rad);
    lon   = rad2deg(lon_rad);
    lat   = rad2deg(lat_rad);
end

end


