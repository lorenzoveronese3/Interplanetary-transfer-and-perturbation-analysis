function plotOrbitalElementsCompact(s)
%PLOTORBITALELEMENTSCOMPACT Compact plot of the evolution of orbital elements.
%
% PROTOTYPE:
%   plotOrbitalElementsCompact(s)
%
% PURPOSE:
%   Creates a single figure (tiled layout) that summarizes over time:
%     a(t), e(t), i(t), RAAN Ω(t) (already unwrapped), argument of perigee ω(t) (already unwrapped).
%
% INPUT:
%   s [struct]  Data to plot (vectors of the same length):
%     - s.t      : time (datetime) [Nx1]
%     - s.a_km   : semi-major axis [km]
%     - s.e      : eccentricity [-]
%     - s.i_deg  : inclination [deg]
%     - s.raan_u : unwrapped RAAN [deg]
%     - s.argp_u : unwrapped argument of perigee [deg]
%     - s.descr  : descriptive string for title/figure name
%
% OUTPUT:
%   None (the function generates a figure).
%
% NOTES:
%   Ω and ω must already be "unwrapped" (unwrap) before calling, to avoid
%   discontinuities due to the 360° modulo.
%
% AUTHORS:
%   Simone Frescura
%   Nicola Pifano
%   Luca Puzzi
%   Lorenzo Veronese
%
%==========================================================================

figure('Name', s.descr);

% Compact 3x2 layout; the last panel (ω) spans two columns for readability.
tl = tiledlayout(3, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
title(tl, s.descr);

% 1) a(t)
nexttile;
plot(s.t, s.a_km, 'LineWidth', 1.2);
grid on; xlabel("Time"); ylabel("a [km]");
title("Semi-major axis a(t)");

% 2) e(t)
nexttile;
plot(s.t, s.e, 'LineWidth', 1.2);
grid on; xlabel("Time"); ylabel("e [-]");
title("Eccentricity e(t)");

% 3) i(t)
nexttile;
plot(s.t, s.i_deg, 'LineWidth', 1.2);
grid on; xlabel("Time"); ylabel("i [deg]");
title("Inclination i(t)");

% 4) Ω(t) (unwrap)
nexttile;
plot(s.t, s.raan_u, 'LineWidth', 1.2);
grid on; xlabel("Time"); ylabel("\Omega [deg] (unwrap)");
title("RAAN \Omega(t) (unwrap)");

% 5) ω(t) (unwrap) across two columns
nexttile([1 2]);
plot(s.t, s.argp_u, 'LineWidth', 1.2);
grid on; xlabel("Time"); ylabel("\omega [deg] (unwrap)");
title("Argument of perigee \omega(t) (unwrap)");

end


