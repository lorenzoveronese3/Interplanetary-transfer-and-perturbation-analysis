function dy = twoBody2BP(~, y, mu)
%TWOBODY2BP ODE system for the two-body problem (2BP) in Cartesian coordinates.
%
% PROTOTYPE:
%   dy = twoBody2BP(t, y, mu)
%
% PURPOSE:
%   Returns the derivative of the Cartesian state y = [r; v] for Keplerian motion
%   about a point-mass central body (2BP model).
%
% DYNAMICAL MODEL:
%   r_dot = v
%   v_dot = -mu * r / ||r||^3
%
% INPUT:
%   t  [1]    Time (not used: autonomous dynamics) [s]
%   y  [6x1]  Cartesian state: [r; v] with r [km], v [km/s]
%   mu [1]    Gravitational parameter of the central body [km^3/s^2]
%
% OUTPUT:
%   dy [6x1]  State derivative: [rdot; vdot] with:
%             rdot = v   [km/s]
%             vdot = a   [km/s^2]
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

% Central gravitational acceleration (Newtonian field, point mass body)
a = -mu * r / (norm(r)^3);  % [km/s^2]

% State derivative for the ODE integrator: [rdot; vdot]
dy = [v; a];
end

