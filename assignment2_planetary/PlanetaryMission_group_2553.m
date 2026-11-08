%% ========================================================================
%  FILE: PlanetaryMission_group_2553.m
%  Planetary Explorer Mission - Assignment 2 (Orbital Mechanics 2025/26)
%  Group: 2553  |  Central body: Venus  |  Assigned perturbations: J2 + SRP
%
% AUTHORS:
%   Simone Frescura
%   Nicola Pifano
%   Luca Puzzi
%   Lorenzo Veronese
%   
%  INPUTS (input data / parameters used):
%   - External files:
%       * 2k_venus_surface.jpg  (RGB image, Venus texture)
%   - Nominal orbit parameters at t0 (hard-coded in BLOCK 1):
%       * a [km], e [-], i [deg], RAAN [deg], omega [deg], theta0 [deg] (theta0=0 as required)
%   - Initial date (needed for SRP) in BLOCK 2c:
%       * t0_str -> conversion to mjd2000_0
%   - Repeating ground track (BLOCK 2b):
%       * k_sat, m_plan (integers)
%   - Perturbed propagation setup (BLOCK 3):
%       * N_orb_pert, Nout_pert, ODE tolerances
%
%  OUTPUTS (results produced by this script):
%   - On-screen figures:
%       * 2BP ground track (BLOCK 2a) for 1T, 3T, 5T
%       * Nominal 3D orbit (1T, 2BP) with scaled Venus (visualization only)
%       * Repeating ground track 2BP (BLOCK 2b)
%       * Ground track with J2+SRP perturbations (BLOCK 2c) for 1T,2T,5T + repeating
%       * Keplerian elements history (BLOCK 4c) and errors (BLOCK 4d)
%       * 3D orbit evolution figure (BLOCK 5)
%
%  DEPENDENCIES (required functions):
%   - Course/WeBeep provided functions:
%       astroConstants, uplanet, date2mjd2000 (in ../course_functions)
%   - Project custom functions (in ./functions):
%       kep2car, car2kep, twoBody2BP, repeatingGroundTrackSemimajorAxis,
%       insertNaNAtDiscontinuity, groundTrackJ2SRP, twoBodyJ2SRP,
%       accPert_J2_SRP_RSW, gauss_eq_motion_rsw, solve_kepler_E
%   - MATLAB built-ins:
%       ode113, odeset, datetime, VideoWriter, exportgraphics
%
% ========================================================================


% ========================================================================
%  GENERAL PURPOSE:
%   Analyze an assigned planetocentric orbit (nominal orbit), study the ground
%   track in both unperturbed and perturbed conditions, propagate the orbit
%   with the assigned perturbations using two different methods and compare
%   them, and finally represent the orbital evolution with an image and a video.
%
%  SCRIPT STRUCTURE (Points 1–5):
%
%   1) NOMINAL ORBIT (BLOCK 1)
%      - Set the nominal Keplerian state at epoch t0 (Initial date).
%      - As required: initial true anomaly theta0 = 0.
%      - Compute main quantities: rp/ra, altitudes hp/ha w.r.t. planet radius,
%        and orbital period T.
%      - Convert the Keplerian state to Cartesian (r0,v0) for the next blocks
%        (ground track and propagations).
%
%   2) GROUND TRACK (BLOCK 2)
%      2a) Nominal orbit in unperturbed 2BP:
%          - Compute the ground track of the nominal trajectory in the two-body
%            model (2BP) for relevant times (multiples of T).
%          - Account for Venus rotation (retrograde) using the convention of the
%            groundTrack function (lon = alpha - thetaG).
%
%      2b) Repeating Ground Track in unperturbed 2BP:
%          - Modify ONLY the semi-major axis a to satisfy the repeating ground
%            track condition with the assigned k:m ratio.
%          - The obtained a_rep is used EXCLUSIVELY for the ground track (as required).
%
%      2c) Ground Track with assigned perturbations:
%          - Repeat the plots of 2a and 2b including the assigned perturbations
%            (J2 + SRP) in the propagation used for the ground track.
%          - The initial epoch ("Initial date") is used for SRP.
%
%   3) ORBIT PROPAGATION WITH PERTURBATIONS (BLOCK 3)
%      IMPORTANT NOTE:
%       - In this block we use the NOMINAL orbit (a from point 1),
%         NOT a_rep (reserved only for ground track in point 2).
%
%      3a) Propagation in Cartesian coordinates (r,v):
%          - Numerical integration of the equations of motion with:
%              a = a_2body + a_J2 + a_SRP
%          - Reconstruct the Keplerian elements history via car2kep
%            to compare against method 3b.
%
%      3b) Propagation in Keplerian elements (Gauss equations, RSW frame):
%          - Numerical integration of Gauss planetary equations.
%          - Perturbing accelerations (J2 + SRP) are projected in the RSW frame
%            and used in Gauss equations.
%          - For Gauss we use the state [a e i RAAN omega M]; the IC given in theta0
%            is converted to M0 (theta0 -> E0 -> M0).
%
%      COMPARISON 3a vs 3b:
%       - Both methods output results at the same time instants using a common
%         time grid t_eval (direct comparison without interpolation).
%
%   4) KEPLERIAN ELEMENTS HISTORY + METHODS COMPARISON (BLOCK 4)
%      - Plot the evolution of a,e,i,Omega,omega,theta versus time.
%      - Plot units: time in days, angles in degrees; unwrap angles to avoid
%        jumps in the plots.
%      - Quantitative comparison: time histories of errors (Gauss - Cartesian) and CPU times.
%
%   5) ORBIT EVOLUTION REPRESENTATION (image + video) (BLOCK 5)
%      - Create a 3D image of the trajectory colored in time, with textured Venus
%        and orbit-start markers (theta wrapping).
%      - Create an MP4 video with accumulated trajectory and current spacecraft position.
%      - This block does NOT recompute dynamics: it uses results from point 3.
%
%  DELIVERY NOTES (from assignment/slides):
%   - initial true anomaly: ALWAYS = 0
%   - "Initial date": epoch t0 associated with the Keplerian state (needed for SRP)
%   - modified semi-major axis a_rep must be used ONLY for ground track
%
%  UNITS USED:
%   - r [km], v [km/s], t [s]  (internally in propagations)
%   - for plots in points 4 and 5: time [days], angles [deg] where specified
%
% If changing the Number of orbit is necessary change the value N_orb_pert
% ========================================================================

clc; clear; close all

% Project functions and data (run this script from its own folder);
% the course-provided functions go in ../course_functions (see README)
addpath('functions', fullfile('..','course_functions'));

%% ========================================================================
%  BLOCK 1 - NOMINAL ORBIT
% ========================================================================

a_km      = 34.1475 * 1e4   ;                        %                     [km]
e         = 0.0142          ;                        %                     [-]
i_deg     = 149             ;                        %                     [deg]
omega_deg = 94              ;                        %                     [deg]
RAAN_deg  = 46              ;                        %                     [deg]
theta_deg = 0               ;                        %                     [deg] 
  
kep0 = [ a_km;
         e;
         deg2rad(i_deg);
         deg2rad(RAAN_deg);
         deg2rad(omega_deg);
         deg2rad(theta_deg) ];

muV = astroConstants(12);                           %                      [km^3/s^2]
RV  = astroConstants(22);                           %                      [km]

rp_km = a_km*(1 - e);                               %                      [km]
ra_km = a_km*(1 + e);                               %                      [km]
hp_km = rp_km - RV;                                 %                      [km]
ha_km = ra_km - RV ;                                %                      [km]
T     = 2*pi*sqrt(a_km^3/muV);                      %                      [s]

[r0, v0] = kep2car(kep0, muV);
y0 = [r0; v0];

%% ========================================================================
%  BLOCK 1 - EXTENSION: Venus rotation and groundTrack convention
% ========================================================================
T_year_V=224.6; % [days]
Trot_s = astroConstants(52) * 3600;   % [s] (retrograde)
omegaP = 2*pi / Trot_s;               % [rad/s] with sign

% groundTrack uses lon = alpha - thetaG.

omegaGT = -omegaP;                   % [rad/s] 

T_h   = T/3600;
T_day = T/86400;

%% ========================================================================
%  PLOTS SETUP (only once) + MAP
% ========================================================================

venusMapFile = fullfile('data','2k_venus_surface.jpg');
venusImg = imread(venusMapFile);

xMapLim  = [-180 180];
yMapLim  = [-90 90];
xMapTick = -180:60:180;
yMapTick = -90:30:90;
mapAlpha = 0.82;

% Default graphic
fsAxes  = 13;
lwAxes  = 1.2;

set(groot, ...
    'defaultFigureColor', 'w', ...
    'defaultAxesFontSize', fsAxes, ...
    'defaultAxesLineWidth', lwAxes, ...
    'defaultAxesXColor', 'k', ...
    'defaultAxesYColor', 'k', ...
    'defaultAxesZColor', 'k', ...
    'defaultAxesBox', 'on', ...
    'defaultAxesTickDir', 'out', ...
    'defaultAxesTickLength', [0.012 0.012], ...
    'defaultLegendBox', 'on', ...
    'defaultLegendColor', 'w', ...
    'defaultLegendEdgeColor', 'k', ...
    'defaultLegendTextColor', 'k', ...
    'defaultLegendFontSize', 11);

lwTrack = 2.0;
mkStart = {'o','LineWidth',2.0,'MarkerSize',8};
mkEnd   = {'s','LineWidth',2.0,'MarkerSize',8};

% Common arguments for groundTrack (2BP)
gtArgs = {'InputType','cart','OutputUnit','deg'};

%% ========================================================================
%  BLOCK 2a - Nominal GROUND TRACK (unperturbed 2BP)
% ========================================================================

tEnd_list   = [1*T, 3*T, 500*T, 9*T];
tLabel_list = {'1T','3T','7T','9T'};
NptsGT = 80000;

for kk = 1:numel(tEnd_list)

    tVec = linspace(0, tEnd_list(kk), NptsGT).';

    [~, ~, lon_deg, lat_deg] = groundTrack(y0, 0, tVec, muV, omegaGT, gtArgs{:});
    [lon_plot, lat_plot] = insertNaNAtDiscontinuity(lon_deg, lat_deg, 180);

    figure;
    hImg = imagesc(xMapLim, yMapLim, venusImg); hImg.AlphaData = mapAlpha;
    set(gca,'YDir','normal','XLim',xMapLim,'YLim',yMapLim,'XTick',xMapTick,'YTick',yMapTick);
    hold on

    if kk == 3
        hTrack = plot(lon_plot, lat_plot, '.', 'MarkerSize', 6);
    else
        hTrack = plot(lon_plot, lat_plot, 'LineWidth', lwTrack);
    end

    hStart = plot(lon_deg(1),   lat_deg(1),  mkStart{:});
    hEnd   = plot(lon_deg(end), lat_deg(end),mkEnd{:});

    xlabel('Longitude [deg]');
    ylabel('Latitude [deg]');
    title(['Nominal ground track - 2BP (Venus) - ', tLabel_list{kk}]);

    legend([hTrack,hStart,hEnd], {'Ground track','Start','End'}, ...
        'Location','southoutside','Orientation','horizontal');
end


%% ========================================================================
%  3D plot of the nominal orbit (only 1T)
%  Note: to make it visible, Venus is scaled ONLY for visualization.
% ========================================================================

optsODE = odeset('RelTol',1e-11,'AbsTol',1e-13);
ode2BP  = @(t,y) twoBody2BP(t,y,muV);

tVec_1T = linspace(0, T, NptsGT).';
[~, ySol_1T] = ode113(ode2BP, tVec_1T, y0, optsODE);

figure('Color','w');

% --- Planet visual scale (ONLY aesthetics) ---
planetScale = 3;     
RV_vis = planetScale*RV;

% Sphere + texture
venusTex = imread(venusMapFile);
[Xe,Ye,Ze] = sphere(80);
Ze = -Ze;  

hSurf = surf(RV_vis*Xe, RV_vis*Ye, RV_vis*Ze);
set(hSurf, ...
    'FaceColor','texturemap', ...
    'CData', venusTex, ...
    'EdgeColor','none', ...
    'FaceLighting','gouraud', ...
    'AmbientStrength', 1);

hold on
plot3(ySol_1T(:,1), ySol_1T(:,2), ySol_1T(:,3), 'LineWidth', 2.5);

axis equal
grid on

xlabel('X [km]');
ylabel('Y [km]');
zlabel('Z [km]');
title('Nominal orbit - 1T');

legend('Venus (scaled)','Nominal orbit (1T)'); 

%% ========================================================================
%  BLOCK 2b - REPEATING GROUND TRACK (unperturbed 2BP)
% ========================================================================

k_sat  = 30;
m_plan = 1;

[a_rep, T_orb_rep, ~] = repeatingGroundTrackSemimajorAxis(k_sat, m_plan, muV, abs(omegaP));

kep_rep = kep0;
kep_rep(1) = a_rep;
kep_rep(6) = 0;

[r0_rep, v0_rep] = kep2car(kep_rep, muV);
y0_rep = [r0_rep; v0_rep];

tEnd_rep = k_sat * T_orb_rep;

NptsRep  = 12000;
tVec_rep = linspace(0, tEnd_rep, NptsRep).';

[~, ~, lon_rep, lat_rep] = groundTrack(y0_rep, 0, tVec_rep, muV, omegaGT, gtArgs{:});
[lon_rep_plot, lat_rep_plot] = insertNaNAtDiscontinuity(lon_rep, lat_rep, 180);

figure;
hImg = imagesc(xMapLim, yMapLim, venusImg); hImg.AlphaData = mapAlpha;
set(gca,'YDir','normal','XLim',xMapLim,'YLim',yMapLim,'XTick',xMapTick,'YTick',yMapTick);
hold on

hTrack = plot(lon_rep_plot, lat_rep_plot, 'LineWidth', lwTrack);
hStart = plot(lon_rep(1),   lat_rep(1),  mkStart{:});
hEnd   = plot(lon_rep(end),lat_rep(end),mkEnd{:});

xlabel('Longitude [deg]');
ylabel('Latitude [deg]');
title('Repeating ground track - 2BP (Venus)');

legend([hTrack,hStart,hEnd], {'Ground track','Start','End'}, ...
    'Location','southoutside','Orientation','horizontal');

%% ========================================================================
%  BLOCK 2c - GROUND TRACK with assigned perturbations: J2 + SRP
% ========================================================================

J2V = astroConstants(32);

Cr       = 1.0;
AoverM   = 5.4766;   % [m^2/kg]
planetID = 2;

t0_str = '2045-8-25-5-24-60';
tmp   = sscanf(t0_str,'%d-%d-%d-%d-%d-%d');
date0 = normalizeDateVector(tmp);     % [Y M D h m s] con carry s->m->h->D
mjd2000_0 = date2mjd2000(date0);


% Common arguments for the perturbed groundTrack
gtPertArgs = { ...
    RV, J2V, mjd2000_0, planetID, Cr, AoverM, ...
    'InputType','cart','OutputUnit','deg','RelTol',1e-11,'AbsTol',1e-13};


% --- (A) Nominal with J2 + SRP ---
tEnd_list_p   = [1*T, 3*T, 7*T, 9*T];
tLabel_list_p = {'1T','3T','7T','9T'};

for kk = 1:numel(tEnd_list_p)

    tVec = linspace(0, tEnd_list_p(kk), NptsGT).';

    [~, ~, lon_deg, lat_deg] = groundTrackJ2SRP(y0, 0, tVec, muV, omegaGT, gtPertArgs{:});
    [lon_plot, lat_plot] = insertNaNAtDiscontinuity(lon_deg, lat_deg, 180);

    figure;
    hImg = imagesc(xMapLim, yMapLim, venusImg); hImg.AlphaData = mapAlpha;
    set(gca,'YDir','normal','XLim',xMapLim,'YLim',yMapLim,'XTick',xMapTick,'YTick',yMapTick);
    hold on

    if kk == 3
        hTrack = plot(lon_plot, lat_plot, '.', 'MarkerSize', 6);
    else
        hTrack = plot(lon_plot, lat_plot, 'LineWidth', lwTrack);
    end

    hStart = plot(lon_deg(1),   lat_deg(1),  mkStart{:});
    hEnd   = plot(lon_deg(end), lat_deg(end),mkEnd{:});

    xlabel('Longitude [deg]');
    ylabel('Latitude [deg]');
    title(['Nominal ground track - J2 + SRP - ', tLabel_list_p{kk}]);

    legend([hTrack,hStart,hEnd], {'Ground track','Start','End'}, ...
        'Location','southoutside','Orientation','horizontal');
end

% --- (B) Repeating with J2 + SRP ---
[~, ~, lon_rep_p, lat_rep_p] = groundTrackJ2SRP(y0_rep, 0, tVec_rep, muV, omegaGT, gtPertArgs{:});
[lon_rep_p_plot, lat_rep_p_plot] = insertNaNAtDiscontinuity(lon_rep_p, lat_rep_p, 180);

figure;
hImg = imagesc(xMapLim, yMapLim, venusImg); hImg.AlphaData = mapAlpha;
set(gca,'YDir','normal','XLim',xMapLim,'YLim',yMapLim,'XTick',xMapTick,'YTick',yMapTick);
hold on

hTrack = plot(lon_rep_p_plot, lat_rep_p_plot, 'LineWidth', lwTrack);
hStart = plot(lon_rep_p(1),   lat_rep_p(1),  mkStart{:});
hEnd   = plot(lon_rep_p(end), lat_rep_p(end),mkEnd{:});

xlabel('Longitude [deg]');
ylabel('Latitude [deg]');
title('Repeating ground track - J2 + SRP (Venus)');

legend([hTrack,hStart,hEnd], {'Ground track','Start','End'}, ...
    'Location','southoutside','Orientation','horizontal');

%% ========================================================================
%  BLOCK 3 - ORBIT PROPAGATION WITH PERTURBATIONS (J2 + SRP)
% ========================================================================

%% --------------------------- 3) Time setup and grid ----------------------
% Practical choice:
%  - propagate for N_orb_pert nominal (Keplerian) orbits.
%  - use a common time grid t_eval to get results directly comparable between
%    3a and 3b (same samples and same instants).
%
% NOTE:
%  - this block uses the NOMINAL orbit (y0 from point 1).
%  - DO NOT use y0_rep (derived from a_rep, only for ground track).

N_orb_pert =500;               % number of orbits to propagate (practical choice)
Tend_pert  = N_orb_pert * T;   % [s] final propagation time

Nout_pert  = 500000;                                   % number of output samples
t_eval     = linspace(0, Tend_pert, Nout_pert).';    % [s] common grid

% Tolerances
optsODE_pert = odeset('RelTol',1e-12,'AbsTol',1e-12);

Nruns_cpu = 10;
doWarmup  = true;

%% ========================================================================
%  BLOCK 3a - Cartesian propagation (r,v)
% ========================================================================
% Integrated state:
%  y = [r; v], with r [km], v [km/s]
%
% Equations of motion:
%  rdot = v
%  vdot = a_2body + a_J2 + a_SRP
%
% Implementation:
%  - external routine twoBodyJ2SRP(t,y,...) returns dy/dt
%  - planet parameters and SRP parameters are those defined in block 2c

odefun_cart_pert = @(t,y) twoBodyJ2SRP( ...
    t, y, muV, RV, J2V, mjd2000_0, planetID, Cr, AoverM);


% Integration on t_eval (ode113 returns the solution exactly at requested samples)
cpu_cart_runs = zeros(Nruns_cpu,1);

if doWarmup
    [~, ~] = ode113(odefun_cart_pert, t_eval, y0, optsODE_pert);
end

for kk = 1:Nruns_cpu
    tic;
    [tC, YC] = ode113(odefun_cart_pert, t_eval, y0, optsODE_pert);
    cpu_cart_runs(kk) = toc;
end

cpu_cart = mean(cpu_cart_runs);
cpu_cart_std = std(cpu_cart_runs);

% In order to get the video faster the number of orbit is set to 100 with
% 5000 points
% N_orb_pert_vid =100;            
% T_end_vid  = N_orb_pert_vid * T;   
% Nout_pert_vid=500;
% 
% t_eval_vid=linspace(0, T_end_vid, Nout_pert_vid).';
% 
% for kk = 1:Nruns_cpu
%     tic;
%     [tC_vid, YC_vid] = ode113(odefun_cart_pert, t_eval, y0, optsODE_pert);
%     cpu_cart_runs(kk) = toc;
% end

% r_vid=YC_vid(:,1:3);   % [km]
% t_days_vid=tC_vid/86400;

% Extract r(t) and v(t)
r_hist = YC(:,1:3);   % [km]
v_hist = YC(:,4:6);   % [km/s]

% Convert to classical Keplerian elements (needed to compare with 3b):
% car2kep returns: [a; e; i; RAAN; omega; theta] with angles in [rad]
kepC_hist = zeros(length(tC), 6);
for k = 1:length(tC)
    kepC_hist(k,:) = car2kep(r_hist(k,:).', v_hist(k,:).', muV).';
end

% Save 3a results in a struct
res3a = struct();
res3a.t        = tC;                 % [s]
res3a.t_days   = tC/86400;           % [days]
res3a.Y        = YC;                 % [km; km/s]
res3a.r        = r_hist;             % [km]
res3a.v        = v_hist;             % [km/s]
res3a.kep      = kepC_hist;          % [a,e,i,RAAN,omega,theta] (angles in rad)
res3a.cpu_time = cpu_cart;           % [s]
res3a.cpu_time_std = cpu_cart_std;   % [s]

%% ========================================================================
%  BLOCK 3b - Propagation with Gauss planetary equations (RSW frame)
% ========================================================================
% Integrated state (standard Gauss formulation):
%  kep = [a; e; i; RAAN; omega; M]
%
% Note:
%  - the project IC is provided in true anomaly theta0 (here theta0 = 0).
%  - Gauss needs the mean anomaly M0:
%      theta0 -> E0 -> M0
%
% Perturbations:
%  - compute total perturbing acceleration (J2 + SRP)
%  - project into RSW frame: [a_R; a_S; a_W]
%  - integrate d/dt of elements with gauss_eq_motion_rsw_safe

% Nominal initial state (from BLOCK 1) for Gauss:
a0  = kep0(1);
e0  = kep0(2);
i0  = kep0(3);
OM0 = kep0(4);
w0  = kep0(5);
th0 = kep0(6);

% Convert theta0 -> E0 -> M0 (elliptic orbit)
E0 = 2*atan( sqrt((1-e0)/(1+e0)) * tan(th0/2) );
E0 = mod(E0, 2*pi);
M0 = E0 - e0*sin(E0);
M0 = mod(M0, 2*pi);

% Initial vector for Gauss
kep0_gauss = [a0; e0; i0; OM0; w0; M0];

% Parameters for the RSW perturbing acceleration routine
pars = struct();
pars.mu           = muV;
pars.R            = RV;
pars.J2           = J2V;
pars.mjd2000_0    = mjd2000_0;
pars.ibody        = planetID;
pars.Cr           = Cr;
pars.AoverM       = AoverM;


% Minimum threshold for e (numerical protection of 1/e terms in Gauss equations)
e_min = 1e-6;

% Total perturbing acceleration (J2 + SRP) in RSW
a_per_rsw = @(t,kep) accPert_J2_SRP_RSW(t, kep, pars);

% Gauss equations (safe version)
odefun_gauss = @(t,kep) gauss_eq_motion_rsw(t, kep, muV, a_per_rsw, e_min);

% Integration on the same time grid t_eval
cpu_gauss_runs = zeros(Nruns_cpu,1);

if doWarmup
    [~, ~] = ode113(odefun_gauss, t_eval, kep0_gauss, optsODE_pert);
end

for kk = 1:Nruns_cpu
    tic;
    [tG, kepG] = ode113(odefun_gauss, t_eval, kep0_gauss, optsODE_pert);
    cpu_gauss_runs(kk) = toc;
end

cpu_gauss = mean(cpu_gauss_runs);
cpu_gauss_std = std(cpu_gauss_runs);

% Reconstruct theta(t) from M(t):
%  - for each time, solve Kepler equation
%  - from E compute theta with a robust formula (atan2)
thetaG = zeros(size(tG));
for k = 1:length(tG)
    e_k = kepG(k,2);
    M_k = kepG(k,6);

    % normalize M in [-pi,pi] (useful for numeric solver)
    Mnorm = mod(M_k + pi, 2*pi) - pi;

    % solve elliptic Kepler equation (external routine)
    E = solve_kepler_E(Mnorm, e_k);

    % E -> theta (robust)
    sinf = (sqrt(1 - e_k^2) * sin(E)) / (1 - e_k*cos(E));
    cosf = (cos(E) - e_k) / (1 - e_k*cos(E));
    th   = atan2(sinf, cosf);
    if th < 0
        th = th + 2*pi;
    end
    thetaG(k) = th;
end

% "Classical" elements for direct comparison with car2kep:
%  kepG_class = [a,e,i,RAAN,omega,theta]
kepG_class = [kepG(:,1), kepG(:,2), kepG(:,3), kepG(:,4), kepG(:,5), thetaG];

% Save 3b results in a struct (useful for point 4 and later)
res3b = struct();
res3b.t         = tG;             % [s]
res3b.t_days    = tG/86400;       % [days]
res3b.kep_gauss = kepG;           % [a,e,i,RAAN,omega,M]
res3b.kep_class = kepG_class;     % [a,e,i,RAAN,omega,theta]
res3b.cpu_time  = cpu_gauss;      % [s]
res3b.cpu_time_std = cpu_gauss_std; % [s]

%% ========================================================================
%  BLOCK 4 - KEPLERIAN ELEMENTS HISTORY + METHODS COMPARISON (Point 4)
% ========================================================================

%% ============================ POINT 4a ================================
Tend_days = res3a.t_days(end);      % total propagation time [days]
Torb_days = Tend_days / N_orb_pert; % "equivalent" orbital period [days]

fprintf('\n==================== POINT 4a ====================\n');
fprintf('Propagation duration: Tend = %.6f days\n', Tend_days);
fprintf('Number of orbits:     N_orb = %d\n', N_orb_pert);
fprintf('Orbital period:       Torb = Tend/N_orb = %.6f days\n', Torb_days);
fprintf('==================================================\n');

%% ============================ POINT 4b ================================
t_days = res3a.t_days;        % [days]

kepC = res3a.kep;             % 3a (Cartesian)
kepG = res3b.kep_class;       % 3b (Gauss)

kepC_u = kepC;
kepG_u = kepG;

kepC_u(:,3) = unwrap(kepC(:,3));
kepC_u(:,4) = unwrap(kepC(:,4));
kepC_u(:,5) = unwrap(kepC(:,5));
kepC_u(:,6) = unwrap(kepC(:,6));

kepG_u(:,3) = unwrap(kepG(:,3));
kepG_u(:,4) = unwrap(kepG(:,4));
kepG_u(:,5) = unwrap(kepG(:,5));
kepG_u(:,6) = unwrap(kepG(:,6));

kepC_deg = kepC_u;
kepG_deg = kepG_u;

kepC_deg(:,3:6) = rad2deg(kepC_u(:,3:6));
kepG_deg(:,3:6) = rad2deg(kepG_u(:,3:6));



%% ============================ POINT 4c ================================
labels = {'a [km]','e [-]','i [deg]','\Omega [deg]','\omega [deg]'};
titles = {'Semi-major axis a(t)', ...
          'Eccentricity e(t)', ...
          'Inclination i(t)', ...
          'RAAN \Omega(t)', ...
          'Argument of periapsis \omega(t)'};

figure('Name','Point 4c - Keplerian elements history','Color','w');

for j = 1:5
    subplot(3,2,j); hold on; grid on;

    plot(t_days, kepC_deg(:,j), 'LineWidth', 1.2);
    plot(t_days, kepG_deg(:,j), '--', 'LineWidth', 1.2);

    xlabel('Time [days]');
    ylabel(labels{j});
    title(titles{j});
end
% Information box
subplot(3,2,6); axis off;
text(0.1, 0.5, ...
    {'\bfElements Evolution:', ...
     '- Blue: Cartesian', ...
     '- Orange: Gauss'}, ...
     'FontSize', 10, 'EdgeColor', 'k', 'BackgroundColor', [1 1 0.9]);
sgtitle('Point 4c - Keplerian elements evolution (J2 + SRP)');

%% ============================ POINT 4d ================================
aC  = kepC_deg(:,1);  eC  = kepC_deg(:,2);  iC  = kepC_deg(:,3);
OC  = kepC_deg(:,4);  wC  = kepC_deg(:,5);  thC = kepC_deg(:,6);

aG  = kepG_deg(:,1);  eG  = kepG_deg(:,2);  iG  = kepG_deg(:,3);
OG  = kepG_deg(:,4);  wG  = kepG_deg(:,5);  thG = kepG_deg(:,6);

da = aG - aC;       % [km]
de = eG - eC;       % [-]
di = iG - iC;       % [deg]

dOmega = mod((OG  - OC ) + 180, 360) - 180;   % [deg]
domega = mod((wG  - wC ) + 180, 360) - 180;   % [deg]
dtheta = mod((thG - thC) + 180, 360) - 180;   % [deg]

errLabels = {'\Delta a [km]','\Delta e [-]','\Delta i [deg]', ...
             '\Delta\Omega [deg]','\Delta\omega [deg]','\Delta\theta [deg]'};

figure('Name','Point 4d - Errors (Gauss - Cartesian)','Color','w');

subplot(3,2,1); grid on;
plot(t_days, da, 'LineWidth', 1.2);
xlabel('Time [days]'); ylabel(errLabels{1}); title('\Delta a(t)');

subplot(3,2,2); grid on;
plot(t_days, de, 'LineWidth', 1.2);
xlabel('Time [days]'); ylabel(errLabels{2}); title('\Delta e(t)');

subplot(3,2,3); grid on;
plot(t_days, di, 'LineWidth', 1.2);
xlabel('Time [days]'); ylabel(errLabels{3}); title('\Delta i(t)');

subplot(3,2,4); grid on;
plot(t_days, dOmega, 'LineWidth', 1.2);
xlabel('Time [days]'); ylabel(errLabels{4}); title('\Delta\Omega(t) (wrap)');

subplot(3,2,5); grid on;
plot(t_days, domega, 'LineWidth', 1.2);
xlabel('Time [days]'); ylabel(errLabels{5}); title('\Delta\omega(t) (wrap)');

% Information box
subplot(3,2,6); axis off;
text(0.1, 0.5, ...
    {'\bfFilter Analysis:', ...
     '- Black: Raw Data', ...
     '- Pink: Long Period', ...
     '- Blue: Secular'}, ...
     'FontSize', 10, 'EdgeColor', 'k', 'BackgroundColor', [1 1 0.9]);


%subplot(3,2,6); grid on;
%plot(t_days, dtheta, 'LineWidth', 1.2);
%xlabel('Time [days]'); ylabel(errLabels{6}); title('\Delta\theta(t) (wrap)');

sgtitle('Point 4d - Errors between methods: Gauss - Cartesian');

rms_local = @(x) sqrt(mean(x.^2));

da_f = da(end);     de_f = de(end);     di_f = di(end);
dO_f = dOmega(end); dw_f = domega(end); dth_f = dtheta(end);

da_max = max(abs(da));     de_max = max(abs(de));     di_max = max(abs(di));
dO_max = max(abs(dOmega)); dw_max = max(abs(domega)); dth_max = max(abs(dtheta));

da_rms = rms_local(da);     de_rms = rms_local(de);     di_rms = rms_local(di);
dO_rms = rms_local(dOmega); dw_rms = rms_local(domega); dth_rms = rms_local(dtheta);

cpu_cart  = res3a.cpu_time;
cpu_gauss = res3b.cpu_time;
cpu_cart_std  = res3a.cpu_time_std;
cpu_gauss_std = res3b.cpu_time_std;

cpu_ratio_runs = cpu_cart_runs ./ cpu_gauss_runs;
cpu_ratio_mean = mean(cpu_ratio_runs);
cpu_ratio_std  = std(cpu_ratio_runs);

fprintf('\n==================== POINT 4d ====================\n');
fprintf('CPU time (mean ± std over %d runs):\n', Nruns_cpu);
fprintf('  3a (Cartesian) = %.6f ± %.6f s\n', cpu_cart, cpu_cart_std);
fprintf('  3b (Gauss)     = %.6f ± %.6f s\n', cpu_gauss, cpu_gauss_std);
fprintf('  Ratio Cart/Gauss = %.3f ± %.3f\n\n', cpu_ratio_mean, cpu_ratio_std);




%% ========================================================================
%  BLOCK 5 - REPRESENT THE EVOLUTION OF THE ORBIT (IMAGE + MOVIE)
% ========================================================================

%% ===================== 5.1) Extract data from res3a ======================
t_days = res3a.t_days(:);   % [days]
r      = res3a.r;           % [km] Nx3
x = r(:,1); y = r(:,2); z = r(:,3);

theta = res3a.kep(:,6);           % [rad]
wrapIdx = find(diff(theta) < -pi);
orbStartIdx = [1; wrapIdx + 1];

%% ===================== 5.2) Venus texture + sphere =======================
tex = imread(venusMapFile);

Rplanet = RV;   % [km]

kPlanet = 5;
Rvis = kPlanet * Rplanet;

nSphere = 180;
[Xs, Ys, Zs] = sphere(nSphere);
Xs = Rvis*Xs; Ys = Rvis*Ys; Zs = Rvis*Zs;

H = size(tex,1); W = size(tex,2);

lon = atan2(Ys, Xs);
lat = asin(Zs / Rvis);
u = (lon + pi) / (2*pi);
v = (lat + pi/2) / pi;

col = max(1, min(W, round(u*(W-1) + 1)));
row = max(1, min(H, round((1-v)*(H-1) + 1)));

CData = zeros(size(Xs,1), size(Xs,2), 3, 'uint8');
for ch = 1:3
    tmp = tex(:,:,ch);
    CData(:,:,ch) = tmp(sub2ind([H W], row, col));
end

%% ===================== 5A) IMAGE: orbit colored in time ==================
figure('Color','w','Name','Point 5 - Orbit evolution (image)');
hold on; grid on; axis equal;

hPlanet_img = surf(Xs, Ys, Zs, CData, ...
    'FaceColor','texturemap', 'EdgeColor','none', 'DisplayName','Venus');


hOrbit_img = scatter3(x, y, z, 2, t_days, 'filled');
hOrbit_img.DisplayName = 'Trajectory (color = time)';

colormap('turbo');
cb = colorbar;
cb.Label.String = 'Time [days]';


xlabel('x [km]'); ylabel('y [km]'); zlabel('z [km]');
title('Point 5 - Orbit evolution (color = time)');
view(3);
rotate3d on;

legend('Location','best');

exportgraphics(gcf, 'point5_orbit_evolution.png', 'Resolution', 300);


%%  6) FILTERING

%   Remove short period oscillation in order to evidence secular evolution
%   and long period parameter variation
% 1. Definition of the filter time window

% Medium acquisition step
dt_days = mean(diff(t_days));

% Acquisition for orbital period
window_size_orb = round(Torb_days / dt_days);
window_size_year=round(T_year_V / dt_days);


% 2. Filter Application

a_raw  = kepG_deg(:,1);
e_raw  = kepG_deg(:,2);
i_raw  = kepG_deg(:,3);
Om_raw = kepG_deg(:,4);
om_raw = kepG_deg(:,5);

a_filt  = movmean(a_raw,  window_size_orb);
a_filt2  = movmean(a_filt,  window_size_year); 
e_filt  = movmean(e_raw,  window_size_orb);
e_filt2  = movmean(e_raw,   window_size_year);
i_filt  = movmean(i_raw,  window_size_orb);
i_filt2  = movmean(i_filt,  window_size_year);
Om_filt = movmean(Om_raw, window_size_orb);
Om_filt2 = movmean(Om_raw,window_size_year);
om_filt = movmean(om_raw, window_size_orb);
om_filt2 = movmean(om_raw, window_size_year);

% 3. Plot (Raw vs Filtered)

figure; grid on; hold on;
plot(t_days, a_raw, 'k', 'LineWidth', 0.5); % Raw (Grigio)
plot(t_days, a_filt, 'm-', 'LineWidth', 1.5); % Filtered (Blu)
plot(t_days, a_filt2, 'b-', 'LineWidth', 1.5);
ylabel('a [km]'); xlabel('Time [days]');
title('Major Semiaxis a(t)'); 
legend
axis tight;
legend('Raw data', 'Long period', 'Secular', 'Location', 'best');



% -- Eccentricity --
figure; hold on; grid on;
plot(t_days, e_raw,  'k', 'LineWidth', 1.5);
plot(t_days, e_filt, 'm-', 'LineWidth', 1.5);
plot(t_days, e_filt2, 'b-', 'LineWidth', 1.5);
ylabel('e [-]'); xlabel('Time [days]');
title('Eccentricity e(t)');
axis tight;
legend('Raw data', 'Long period', 'Secular', 'Location', 'best');

% inclination
figure; hold on; grid on;
plot(t_days, i_raw,  'k', 'LineWidth', 0.5);
plot(t_days, i_filt, 'm-', 'LineWidth', 1.5);
plot(t_days, i_filt2, 'b-', 'LineWidth', 1.5);
ylabel('i [deg]'); xlabel('Time [days]');
title('Inclination i(t)');
axis tight;
legend('Raw data', 'Long period', 'Secular', 'Location', 'best');

% RAAN
figure; hold on; grid on;
plot(t_days, Om_raw, 'k',  'LineWidth', 0.5);
plot(t_days, Om_filt,'m-', 'LineWidth', 1.5);
plot(t_days, Om_filt2,'b-', 'LineWidth', 1.5);
ylabel('\Omega [deg]'); xlabel('Time [days]');
title('RAAN \Omega(t) (Unwrapped)');
axis tight;
legend('Raw data', 'Long period', 'Secular', 'Location', 'best');

% Argument of pericentre
figure; hold on; grid on;
plot(t_days, om_raw, 'k', 'LineWidth', 0.5);
plot(t_days, om_filt,'m-', 'LineWidth', 1.5);
plot(t_days, om_filt2,'b-', 'LineWidth', 1.5);
ylabel('\omega [deg]'); xlabel('Time [days]');
title('Arg. of Pericenter \omega(t) (Unwrapped)');
axis tight;
legend('Raw data', 'Long period', 'Secular', 'Location', 'best');

% Definiamo i dati
x = 0:0.1:10;
y1 = sin(x);
y2 = x.^2;

 



%% Box filtered Plot 


figure('Name', 'Filtered vs Raw Evolution', 'Color', 'w');

% -- Major Semiaxis --
subplot(3,2,1); hold on; grid on;
plot(t_days, a_raw,  'k', 'LineWidth', 0.5); 
plot(t_days, a_filt, 'm-', 'LineWidth', 1.5); 
plot(t_days, a_filt2, 'b-', 'LineWidth', 1.5); 
ylabel('a [km]'); xlabel('Time [days]');
title('Major Semiaxis a(t)');
axis tight;

% -- Eccenticity --
subplot(3,2,2); hold on; grid on;
plot(t_days, e_raw,  'k',  'LineWidth', 0.5);
plot(t_days, e_filt, 'm-', 'LineWidth', 1.5);
plot(t_days, e_filt2, 'b-', 'LineWidth', 1.5);
ylabel('e [-]'); xlabel('Time [days]');
title('Eccentricity e(t)');
axis tight;

% -- Inclination --
subplot(3,2,3); hold on; grid on;
plot(t_days, i_raw,  'k', 'LineWidth', 0.5);
plot(t_days, i_filt, 'm-', 'LineWidth', 1.5);
plot(t_days, i_filt2, 'b-', 'LineWidth', 1.5);
ylabel('i [deg]'); xlabel('Time [days]');
title('Inclination i(t)');
axis tight;

% -- RAAN (Omega) --
subplot(3,2,4); hold on; grid on;
plot(t_days, Om_raw, 'k', 'LineWidth', 0.5);
plot(t_days, Om_filt,'m-', 'LineWidth', 1.5);
plot(t_days, Om_filt2,'b-', 'LineWidth', 1.5);
ylabel('\Omega [deg]'); xlabel('Time [days]');
title('RAAN \Omega(t) (Unwrapped)');
axis tight;

% -- Argomento Pericentro (omega) --
subplot(3,2,5); hold on; grid on;
plot(t_days, om_raw, 'k', 'LineWidth', 0.5);
plot(t_days, om_filt,'m-', 'LineWidth', 1.5);
plot(t_days, om_filt2,'b-', 'LineWidth', 1.5);
ylabel('\omega [deg]'); xlabel('Time [days]');
title('Arg. Pericentro \omega(t) (Unwrapped)');
axis tight;

% Information box
subplot(3,2,6); axis off;
text(0.1, 0.5, ...
    {'\bfFilter Analysis:', ...
     '- Black: Raw Data', ...
     '- Pink: Long Period', ...
     '- Blue: Secular'}, ...
     'FontSize', 10, 'EdgeColor', 'k', 'BackgroundColor', [1 1 0.9]);

% sgtitle(['Punto 6: Filtraggio Alte Frequenze (Finestra \approx 1 Periodo Orbitale)']);
% 
% fprintf('Grafici filtrati generati.\n');
% fprintf('==================================================\n');


%% ========================================================================
%  FILE: RealDataAnalysis.m
%  Point 7 - Real data analysis (LEO + MEO) 
%
%  PURPOSE:
%   Analyze real data (historical TLEs) for two objects in different orbital
%   regions (LEO and MEO) and:
%     1) build ordered time series of orbital elements,
%     2) plot the evolution of a, e, i, Ω, ω,
%     3) estimate the mean precession rates of Ω and ω via linear fitting.
%
%  ANALYZED OBJECTS:
%   - LEO : NORAD 13819
%   - MEO : NORAD 8820 (LAGEOS-1)
%
%  INPUT (files in the working folder):
%   - Historical TLE file LEO: "LEO_13819*" (text with TLE line pairs)
%   - Historical TLE file MEO: "MEO_8820*" (text with TLE line pairs)
%
%  OUTPUT:
%   - Console prints:
%       * number of TLEs, time span, mean rates dΩ/dt and dω/dt [deg/day]
%   - Figures:
%       * one figure per object with a(t), e(t), i(t), Ω(t) unwrap, ω(t) unwrap
%
%  DEPENDENCIES:
%   - readHistoricalTLE.m
%   - plotOrbitalElementsCompact.m
%   - MATLAB: dir, sort, polyfit, unwrap, deg2rad, rad2deg
%
%  NOTES:
%   - TLEs contain "mean" elements (consistent with the SGP model).
%     Here we perform a descriptive trend analysis, not numerical propagation.
%   - Unwrap of Ω and ω is required to remove discontinuities due to the 360° modulo.
%
% ========================================================================

%clear; clc; close all;

%% ------------------------------------------------------------------------
% 1) Object definition and file patterns
%    (it is assumed that the folder contains at least one file for each pattern)
% -------------------------------------------------------------------------
obj(1).tag      = "LEO";
obj(1).norad    = 13819;
obj(1).pattern  = "LEO_13819*";
obj(1).descr    = "LEO (13819) - Object in low Earth orbit";

obj(2).tag      = "MEO";
obj(2).norad    = 8820;
obj(2).pattern  = "MEO_8820*";
obj(2).descr    = "MEO (8820) - LAGEOS-1";

%% ------------------------------------------------------------------------
% 2) Read TLE + build ordered time series + linear fit
% -------------------------------------------------------------------------
for k = 1:numel(obj)

    % File selection: first match of the pattern 
    files = dir(fullfile('data', obj(k).pattern));
    if isempty(files)
        error('TLE file %s not found in the data folder: download it from Space-Track (see README).', obj(k).pattern);
    end
    filePath = fullfile(files(1).folder, files(1).name);
    obj(k).file = filePath;  

    % Read TLE (parsing + computation of a, altitudes, etc.)
    tle = readHistoricalTLE(filePath);

    % Sort samples by time (required for consistent plots and fits)
    [t_sorted, idx] = sort(tle.epoch);

    obj(k).t        = t_sorted;
    obj(k).a_km     = tle.a_km(idx);
    obj(k).e        = tle.e(idx);
    obj(k).i_deg    = tle.i_deg(idx);
    obj(k).raan_deg = tle.raan_deg(idx);
    obj(k).argp_deg = tle.argp_deg(idx);

    % Numeric time axis for the fit: days from the initial instant
    obj(k).t0      = obj(k).t(1);
    obj(k).dt_days = days(obj(k).t - obj(k).t0);

    % ---------------------------------------------------------------------
    % Angle unwrap (in degrees) to remove 360 -> 0 jumps:
    %   degrees -> radians -> unwrap -> degrees
    % ---------------------------------------------------------------------
    obj(k).raan_u = rad2deg(unwrap(deg2rad(obj(k).raan_deg)));
    obj(k).argp_u = rad2deg(unwrap(deg2rad(obj(k).argp_deg)));

    % Linear fit: estimate mean rates [deg/day]
    p_raan = polyfit(obj(k).dt_days, obj(k).raan_u, 1);
    p_argp = polyfit(obj(k).dt_days, obj(k).argp_u, 1);
    obj(k).dRAAN_dt = p_raan(1);
    obj(k).dARGP_dt = p_argp(1);

    % ---------------------------------------------------------------------
    % Summary print (without "File used: ...", as requested)
    % ---------------------------------------------------------------------
    fprintf("\n============================================================\n");
    fprintf("%s\n", obj(k).descr);
    fprintf("Number of TLEs: %d\n", numel(obj(k).t));
    fprintf("Time span: %s  ->  %s\n", string(obj(k).t(1)), string(obj(k).t(end)));
    fprintf("Estimated mean rates (linear fit):\n");
    fprintf("  dΩ/dt  ≈ %+.6f deg/day\n", obj(k).dRAAN_dt);
    fprintf("  dω/dt  ≈ %+.6f deg/day\n", obj(k).dARGP_dt);
    fprintf("============================================================\n");
end

%% ------------------------------------------------------------------------
% 3) Plots: one "compact" figure for each object
% -------------------------------------------------------------------------
for k = 1:numel(obj)
    plotOrbitalElementsCompact(obj(k));
end


