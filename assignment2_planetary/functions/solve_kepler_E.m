function E = solve_kepler_E(M, e)
%SOLVE_KEPLER_E Solves the elliptical Kepler equation: M = E - e*sin(E).
%
% PROTOTYPE:
%   E = solve_kepler_E(M, e)
%
% PURPOSE:
%   Given mean anomaly M and eccentricity e (elliptical orbit), computes
%   the eccentric anomaly E via an iterative Newton-Raphson method.
%
% INPUT:
%   M [1]  Mean anomaly [rad] (typically normalized in [-pi, pi] for numerical stability)
%   e [1]  Eccentricity [-]
%
% OUTPUT:
%   E [1]  Eccentric anomaly [rad]
%
% METHOD:
%   Newton-Raphson applied to the function:
%     f(E)  = E - e*sin(E) - M
%     f'(E) = 1 - e*cos(E)
%   with a simple but robust initial guess.
%
% NOTES:
%   - The stopping criterion is based on the increment |dE|.
%   - If convergence is not achieved within maxIter, the last estimate is returned.
%
% AUTHORS:
%   Simone Frescura
%   Nicola Pifano
%   Luca Puzzi
%   Lorenzo Veronese
%==========================================================================

% --------------------- Initial guess (robust) ----------------------------
% For small eccentricities E ≈ M; for high eccentricities a guess
% farther from M can help Newton converge.
if e < 0.8
    E = M;
else
    E = pi * sign(M);
    if E == 0
        E = pi;
    end
end

% --------------------- Iteration parameters ------------------------------
maxIter = 50;
tol     = 1e-12;

% --------------------- Newton-Raphson -----------------------------------
for k = 1:maxIter
    f  = E - e*sin(E) - M;
    fp = 1 - e*cos(E);

    dE = -f/fp;
    E  = E + dE;

    if abs(dE) < tol
        break;
    end
end

end

