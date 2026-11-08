function kep = car2kep(r, v, mu)

%CAR2KEP Converts Cartesian state into Keplerian elements.
%
%   kep = car2kep(r, v, mu)
%
%   Input:
%       r  : 3×1 position vector in the inertial frame [km]
%       v  : 3×1 velocity vector in the inertial frame [km/s]
%       mu : gravitational parameter of the central body [km^3/s^2]
%
%   Output:
%       kep : 6×1 vector of Keplerian elements
%             kep = [a; e; i; RAAN; omega; theta];
%             a     = semi-major axis [km]
%                     (WARNING: in the parabolic case e≈1, here a contains p)
%             e     = eccentricity [-]
%             i     = inclination [rad]
%             RAAN  = Ω, right ascension of the ascending node [rad]
%             omega = argument of periapsis [rad]
%             theta = true anomaly [rad]
%
%   The function handles:
%       - elliptical orbits (0 ≤ e < 1)
%       - nearly circular orbits (very small e)
%       - parabolic orbits (e ≈ 1)
%       - hyperbolic orbits (e > 1)
%
%   Conventions and special cases:
%       - Inertial IJK frame: the same used in kep2car.
%       - Perifocal frame (p-q-w): (p-q) plane = orbital plane, p toward
%         periapsis, q rotated by +90° in the plane, w perpendicular to the plane.
%       - For e ≈ 1 (parabolic) it is assumed that the orbital parameter
%             p = h^2 / mu
%         is returned in kep(1) instead of the semi-major axis (which
%         would be infinite). This is consistent with the convention used
%         in the kep2car function.
%       - For nearly circular (e < 1e-10) and/or equatorial (i ≈ 0) orbits
%         RAAN and omega become ill-defined. The following choices are used:
%           * e < tol, i > tol  → inclined circular orbit:
%               RAAN from node vector,
%               omega = 0,
%               theta = argument of latitude u.
%           * e < tol, i ≈ 0   → equatorial circular orbit:
%               RAAN = 0,
%               omega = 0,
%               theta = true longitude (angle in the equatorial plane).
%           * e > tol, i ≈ 0   → equatorial non-circular orbit:
%               RAAN = 0,
%               omega = longitude of periapsis.
%
%   Numerical notes:
%       - All angles are reported in the interval [0, 2π).
%       - Tolerances are used to distinguish degenerate cases:
%           e_circ_tol  ~ 1e-10
%           i_eq_tol    ~ 1e-10
%           e_parab_tol ~ 1e-8
%
% AUTHORS:
%   Simone Frescura
%   Nicola Pifano
%   Luca Puzzi
%   Lorenzo Veronese
%==========================================================================
   

    r = r(:);
    v = v(:);

   

    r_norm = norm(r);
    v_norm = norm(v);

    

    %% Numerical tolerances
    e_circ_tol  = 1e-10;  % threshold for nearly circular orbit
    i_eq_tol    = 1e-10;  % threshold for nearly equatorial orbit (i~0)
    e_parab_tol = 1e-8;   % tolerance for e ≈ 1 (parabolic)
    n_tol       = 1e-12;  % threshold for null node vector

    k_hat = [0; 0; 1];

    %% Fundamental vector quantities
    % Specific angular momentum
    h = cross(r, v);
    h_norm = norm(h);
    
    h_hat = h / h_norm;

    % Ascending node vector: n = k × h
    n = cross(k_hat, h);
    n_norm = norm(n);

    % Eccentricity vector
    e_vec = (cross(v, h) / mu) - (r / r_norm);
    e = norm(e_vec);

    % Inclination
    i = acos( h(3) / h_norm );

    %% Semi-major axis (or parameter p for parabolic orbits)
    % Specific energy
    eps = 0.5 * v_norm^2 - mu / r_norm;

    if abs(e - 1.0) <= e_parab_tol
        % Parabolic orbit: eps ≈ 0, a → ∞.
        % Use p = h^2 / mu and return a = p for compatibility with kep2car.
        p = h_norm^2 / mu;
        a = p;
        e = 1.0;  % force e = 1 for numerical consistency
    else
        % Elliptical, circular, hyperbolic orbits (e ≠ 1)
        a = -mu / (2 * eps);
    end

    %% RAAN angle (right ascension of the ascending node)
    if n_norm > n_tol
        n_hat = n / n_norm;
        RAAN = atan2(n_hat(2), n_hat(1));
        if RAAN < 0
            RAAN = RAAN + 2*pi;
        end
    else
        % Equatorial orbit (i ≈ 0): ascending node not defined
        RAAN = 0;
    end

    %% Argument of periapsis, true anomaly, and degenerate cases

    % Useful unit directions
    r_hat = r / r_norm;

    if e > e_circ_tol && n_norm > n_tol && i > i_eq_tol
        % General case: non-circular and non-equatorial orbit

        e_hat = e_vec / e;
        n_hat = n / n_norm;

        % Argument of periapsis omega
        cos_omega = dot(n_hat, e_hat);
        sin_omega = dot(cross(n_hat, e_hat), h_hat);
        omega = atan2(sin_omega, cos_omega);
        if omega < 0
            omega = omega + 2*pi;
        end

        % True anomaly theta
        cos_theta = dot(e_hat, r_hat);
        sin_theta = dot(cross(e_hat, r_hat), h_hat);
        theta = atan2(sin_theta, cos_theta);
        if theta < 0
            theta = theta + 2*pi;
        end

    elseif e <= e_circ_tol && n_norm > n_tol && i > i_eq_tol
        % Nearly circular inclined orbit:
        % periapsis is not defined, we use:
        %   omega = 0
        %   theta = argument of latitude u (measured from ascending node)

        e = 0;  % force e = 0 for consistency

        n_hat = n / n_norm;

        cos_u = dot(n_hat, r_hat);
        sin_u = dot(cross(n_hat, r_hat), h_hat);
        u = atan2(sin_u, cos_u);
        if u < 0
            u = u + 2*pi;
        end

        omega = 0;
        theta = u;

    elseif e > e_circ_tol && (i <= i_eq_tol || n_norm <= n_tol)
        % Equatorial non-circular orbit:
        % ascending node undefined, use:
        %   RAAN  = 0
        %   omega = longitude of periapsis (varpi)
        %   theta = true anomaly in the equatorial plane

        RAAN = 0;

        e_hat = e_vec / e;

        % Longitude of periapsis (angle from inertial x-axis to e_vec)
        varpi = atan2(e_hat(2), e_hat(1));
        if varpi < 0
            varpi = varpi + 2*pi;
        end
        omega = varpi;

        % True anomaly as the angle between e_vec and r in the equatorial plane
        cos_theta = dot(e_hat, r_hat);
        sin_theta = dot(cross(e_hat, r_hat), h_hat);
        theta = atan2(sin_theta, cos_theta);
        if theta < 0
            theta = theta + 2*pi;
        end

    else
        % Nearly circular equatorial orbit:
        %   e ~ 0, i ~ 0 → RAAN and omega not defined.
        % Use:
        %   RAAN  = 0
        %   omega = 0
        %   theta = true longitude (angle between x-axis and r in the equatorial plane)

        e = 0;
        RAAN = 0;
        omega = 0;

        theta = atan2(r_hat(2), r_hat(1));
        if theta < 0
            theta = theta + 2*pi;
        end
    end

    %% Output
    kep = [a; e; i; RAAN; omega; theta];
end


