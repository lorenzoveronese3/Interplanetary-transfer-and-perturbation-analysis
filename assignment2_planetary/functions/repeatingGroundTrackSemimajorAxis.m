function [a, T, n] = repeatingGroundTrackSemimajorAxis(k, m, mu, omegaE)
%REPEATINGGROUNDTRACKSEMIMAJORAXIS Semi-major axis for a repeating ground track (2BP model).
%
% PROTOTYPE:
%   [a, T, n] = repeatingGroundTrackSemimajorAxis(k, m, mu, omegaE)
%
% PURPOSE:
%   Computes the semi-major axis a (and the corresponding period T) that satisfies the
%   repeating ground track condition in the ideal model:
%     - Keplerian orbit (2BP)
%     - central body rotating at constant angular velocity omegaE
%
% REPEAT CONDITION (k : m):
%   After k orbital revolutions of the satellite, the planet completes m full rotations.
%   In terms of periods:
%       k * T = m * Trot
%   with:
%       Trot = 2*pi / omegaE
%   therefore:
%       T = (m/k) * (2*pi/omegaE)
%
% KEPLER'S LAW (2BP):
%       T = 2*pi * sqrt(a^3/mu)
%   from which:
%       a = ( mu * (T/(2*pi))^2 )^(1/3)
%
% INPUT:
%   k      [1]   Number of orbits after which the ground track repeats [-]
%   m      [1]   Number of planet rotations in the same interval [-]
%   mu     [1]   Gravitational parameter of the central body [km^3/s^2]
%   omegaE [1]   Planet rotation angular velocity [rad/s]
%
% OUTPUT:
%   a      [1]   Required semi-major axis [km]
%   T      [1]   Required orbital period [s]
%   n      [1]   Mean motion [rad/s], n = 2*pi/T
%
% NOTES ON UNITS:
%   The units of a follow those of mu (here: km^3/s^2 -> a in km). Time is in seconds.
%
% AUTHORS:
%   Simone Frescura
%   Nicola Pifano
%   Luca Puzzi
%   Lorenzo Veronese
%==========================================================================


% Planet rotation period (uniform rotation)
Trot = 2*pi / omegaE;     % [s]

% Orbital period satisfying the repetition ratio k:m
T = Trot * (m / k);       % [s]

% Semi-major axis from Kepler's third law (2BP)
a = ( mu * (T/(2*pi))^2 )^(1/3);   % [km]

% Mean motion
n = 2*pi / T;             % [rad/s]

end

