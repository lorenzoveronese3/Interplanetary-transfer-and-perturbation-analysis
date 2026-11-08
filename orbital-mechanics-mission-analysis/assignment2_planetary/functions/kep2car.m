function [r, v] = kep2car(kep, mu)


%KEP2CAR Converts Keplerian elements into Cartesian state.
%
%   [r, v] = kep2car(kep, mu)
%
%   Input:
%       kep : 6×1 column vector with Keplerian elements
%             kep = [a; e; i; RAAN; omega; theta];
%             a     = semi-major axis [km]
%             e     = eccentricity [-]
%             i     = inclination [rad]
%             RAAN  = Ω, right ascension of the ascending node [rad]
%             omega = argument of periapsis [rad]
%             theta = true anomaly [rad]
%
%       mu  : gravitational parameter of the central body [km^3/s^2]
%
%   Output:
%       r : 3×1 column vector, Cartesian position in the inertial frame [km]
%       v : 3×1 column vector, Cartesian velocity in the inertial frame [km/s]
%
%   The function handles:
%       - elliptical orbits (0 ≤ e < 1)
%       - nearly circular orbits (very small e)
%       - parabolic orbits (e ≈ 1)
%       - hyperbolic orbits (e > 1)
%
%   Conventions:
%       - Perifocal frame (p-q-w): the (p-q) plane coincides with the orbital plane,
%         the p axis points toward periapsis, q is rotated by +90° in the orbital plane,
%         w is perpendicular to the orbital plane.
%       - Inertial frame (I-J-K): inertial reference frame (e.g., ECI for Earth).
%       - The transformation from perifocal to inertial is given by:
%           Q = R3(RAAN) * R1(i) * R3(omega)
%         where R3 is a rotation about the z-axis, R1 a rotation about the x-axis.
%
%   Note on special cases:
%       - If e is very small (e < 1e-10), it is treated as 0 (circular orbit)
%         to improve numerical stability.
%       - If e ≈ 1 (within a small tolerance), a parabolic orbit is assumed.
%         In this case it is ASSUMED that the first element 'a' of kep represents
%         directly the orbital parameter p (semi-latus rectum) instead of the semi-major axis:
%             p = kep(1)
%         This avoids issues with a → ∞ for parabolic orbits.
%       - For hyperbolic orbits (e > 1) the same position/velocity formula is used
%         with p = a(1 - e^2); for the speed magnitude we use sqrt(mu/abs(p))
%         to handle possible inconsistent signs of p (e.g., non-standard sign of a).
%
% AUTHORS:
%   Simone Frescura
%   Nicola Pifano
%   Luca Puzzi
%   Lorenzo Veronese
%%-----------------------------------------------------------------------------------------

    

    % Force kep to be a column (6x1)
    kep = kep(:);

   
    % Extract Keplerian elements
    a     = kep(1);  % semi-major axis [km] (or p in the parabolic case, see below)
    e     = kep(2);  % eccentricity
    inc   = kep(3);  % inclination [rad]
    RAAN  = kep(4);  % right ascension of the ascending node Ω [rad]
    omega = kep(5);  % argument of periapsis [rad]
    theta = kep(6);  % true anomaly [rad]

    % Numerical tolerances
    e_circ_tol = 1e-10;   % threshold to consider the orbit circular
    e_parab_tol = 1e-8;   % tolerance to consider e ≈ 1 (parabolic)

    % "Effective" copy of e for possible numerical corrections
    e_eff = e;

    %% Orbit classification and computation of the orbital parameter p
    %
    % p is the orbital parameter (semi-latus rectum):
    %   - for elliptical/circular orbits: p = a(1 - e^2)
    %   - for hyperbolic orbits: same (note that usually a < 0)
    %   - for parabolic orbits: a is formally infinite, so p is used directly
    %     provided by the user instead of a.

    if abs(e) < e_circ_tol
        % Nearly circular orbit → treated as circular
        e_eff = 0.0;
        p = a * (1 - e_eff^2);

        if p <= 0
            error('kep2car:ParamError', ...
                ['Non-positive orbital parameter p in the nearly circular case. ', ...
                 'Check a and mu.']);
        end

    elseif abs(e - 1.0) <= e_parab_tol
        % Parabolic orbit (e ≈ 1)
        e_eff = 1.0;

        % Documented assumption: the user provides p in the first element of kep.
        % Therefore:
        %   p = kep(1);
        % and not the semi-major axis (which would be infinite for e = 1).
        p = a;

        if p <= 0
            error('kep2car:ParabolicParamError', ...
                ['For parabolic orbits it is assumed that the first element of kep ', ...
                 '(a) represents the orbital parameter p = a > 0. ', ...
                 'Received value: p = %g.'], p);
        end

    elseif e < 1.0
        % Generic elliptical orbit (0 < e < 1)
        p = a * (1 - e_eff^2);

        if p <= 0
            error('kep2car:EllipseParamError', ...
                ['Non-positive orbital parameter p = a(1-e^2). ', ...
                 'Check the sign and value of a and e.']);
        end

    else
        % Hyperbolic orbit (e > 1)
        p = a * (1 - e_eff^2);
        % For hyperbolic orbits the user should use a < 0, so p is positive.
        % However, for NUMERICAL robustness in velocity we will use abs(p) in the denominator
        % of sqrt(mu/abs(p)), as specified.

        if p == 0
            error('kep2car:HyperbolaParamError', ...
                ['Orbital parameter p = a(1-e^2) is zero for a hyperbolic orbit. ', ...
                 'Check the input.']);
        end
    end

    %% Position and velocity vectors in the perifocal frame
    %
    % Common formula (elliptical, circular, parabolic, hyperbolic):
    %
    % r_pf = [ p*cos(theta)/(1 + e*cos(theta));
    %          p*sin(theta)/(1 + e*cos(theta));
    %          0 ];
    %
    % v_pf = sqrt(mu/p) * [ -sin(theta);
    %                        e + cos(theta);
    %                        0 ];
    %
    % For hyperbolic (e>1) we use sqrt(mu/abs(p)) to avoid issues
    % if p turned out negative due to "non-standard" choices of a.

    denom = 1 + e_eff * cos(theta);
    

    r_pf = [ p * cos(theta) / denom;
             p * sin(theta) / denom;
             0 ];

    if e_eff > 1.0
        % Hyperbolic: use abs(p) in the velocity
        v_pf = sqrt(mu / abs(p)) * [ -sin(theta);
                                      e_eff + cos(theta);
                                      0 ];
    else
        % Elliptical / Circular / Parabolic
        v_pf = sqrt(mu / p) * [ -sin(theta);
                                 e_eff + cos(theta);
                                 0 ];
    end

    %% Rotation matrix from perifocal frame to inertial frame
    %
    % The rotations are:
    %   1) Rotation R3_Om about the z-axis by RAAN (Ω)
    %   2) Rotation R1_i  about the x-axis by inclination i
    %   3) Rotation R3_om about the z-axis by omega (ω)
    %
    % Q = R3_Om * R1_i * R3_om
    %
    % Thus:
    %   [r_IJK] = Q * [r_pqw]
    %   [v_IJK] = Q * [v_pqw]

    R3_Om = [ cos(RAAN)  -sin(RAAN)  0;
              sin(RAAN)   cos(RAAN)  0;
              0           0          1 ];

    R1_i = [ 1      0           0;
             0    cos(inc)   -sin(inc);
             0    sin(inc)    cos(inc) ];

    R3_om = [ cos(omega)  -sin(omega)  0;
              sin(omega)   cos(omega)  0;
              0            0           1 ];

    Q = R3_Om * R1_i * R3_om;

    %% Position and velocity in the inertial frame
    r = Q * r_pf;
    v = Q * v_pf;
end



