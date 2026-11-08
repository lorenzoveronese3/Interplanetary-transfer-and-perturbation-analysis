function tle = readHistoricalTLE(filename)
%READHISTORICALTLE Reading historical TLEs (2-line format) and reconstruction of orbital parameters.
%
% PROTOTYPE:
%   tle = readHistoricalTLE(filename)
%
% PURPOSE:
%   Reads a text file containing multiple TLE pairs (line 1 and line 2) and
%   returns, for each epoch, the main elements contained in the TLE.
%   From the mean motion n (rev/day) it also derives the semi-major axis a via
%   Kepler's 3rd law (2BP model).
%
% INPUT:
%   filename [1]  Path/name of the text file with TLEs in standard format
%
% OUTPUT (struct tle):
%   tle.norad_id   [1]    NORAD ID (integer)
%   tle.epoch      [Nx1]  TLE epoch (datetime)
%   tle.i_deg      [Nx1]  Inclination [deg]
%   tle.raan_deg   [Nx1]  RAAN (Omega) [deg]
%   tle.e          [Nx1]  Eccentricity [-]
%   tle.argp_deg   [Nx1]  Argument of perigee (omega) [deg]
%   tle.M_deg      [Nx1]  Mean anomaly [deg]
%   tle.n_rev_day  [Nx1]  Mean motion [rev/day]
%   tle.a_km       [Nx1]  Semi-major axis [km]
%   tle.h_mean_km  [Nx1]  Mean altitude (a - R_E) [km]
%   tle.hp_km      [Nx1]  Perigee altitude [km]
%   tle.ha_km      [Nx1]  Apogee altitude [km]
%
% NOTES:
%   - The file is assumed to be consistent: alternating TLE lines "1 " and "2 " .
%   - The TLE epoch is in the format YYDDD.DDDDDDDD (day of year with fraction).
%   - Constants: standard Earth gravitational parameter and Earth equatorial radius.
%
% AUTHORS:
%   - Simone Frescura
%   - Nicola Pifano
%   - Luca Puzzi
%   - Lorenzo Veronese 
%
%==========================================================================

% Constants (consistent units: km, s)
mu = 398600.4418; % [km^3/s^2]
Re = 6378.137;    % [km]

% ---------------------------- File reading -------------------------------
txt   = fileread(char(filename));
lines = splitlines(string(txt));
lines = strtrim(lines);
lines(lines == "") = [];

% --------------------- Extract TLE pairs (1/2) ---------------------------
L1 = strings(0,1);
L2 = strings(0,1);

tempL1 = "";
for k = 1:numel(lines)
    r = lines(k);

    if startsWith(r, "1 ")
        tempL1 = r;

    elseif startsWith(r, "2 ")
        % Assume consistent file: each "2 " has the previous "1 "
        L1(end+1,1) = tempL1; 
        L2(end+1,1) = r;      
        tempL1 = "";
    end
end

N = numel(L1);

% ----------------------- Output preallocation ----------------------------
epoch  = NaT(N,1);
i_deg  = zeros(N,1);
raan   = zeros(N,1);
ecc    = zeros(N,1);
argp   = zeros(N,1);
Mdeg   = zeros(N,1);
nrev   = zeros(N,1);

% NORAD ID from the first line 1 (columns 3–7)
norad_id = str2double(extractBetween(L1(1), 3, 7));

% -------------------------- Fixed-column parsing -------------------------
for j = 1:N
    a1 = char(L1(j));
    a2 = char(L2(j));

    % ---- Epoch: YYDDD.DDDDDDDD (line 1, columns 19–32)
    % YY = last two digits of the year; DDD.D... = day of year with fraction
    epNum = str2double(strtrim(string(a1(19:32))));
    yy  = floor(epNum/1000);
    doy = epNum - yy*1000;

    if yy < 57
        Y = 2000 + yy;
    else
        Y = 1900 + yy;
    end

    % DOY -> datetime conversion: reference to Jan 1 + (doy-1) days
    epoch(j) = datetime(Y,1,1,0,0,0) + days(doy - 1);

    % ---- Orbital elements (line 2, standard columns)
    i_deg(j) = str2double(strtrim(string(a2(9:16))));
    raan(j)  = str2double(strtrim(string(a2(18:25))));

    eccStr = strtrim(string(a2(27:33))); % without "0."
    ecc(j) = str2double("0." + eccStr);

    argp(j) = str2double(strtrim(string(a2(35:42))));
    Mdeg(j) = str2double(strtrim(string(a2(44:51))));
    nrev(j) = str2double(strtrim(string(a2(53:63))));
end

% --------------------- Kepler: from n to semi-major axis -----------------
% n [rev/day] -> [rad/s], then a = (mu/n^2)^(1/3)
n_rad_s = nrev * 2*pi / 86400;          % [rad/s]
a_km    = (mu ./ (n_rad_s.^2)).^(1/3);  % [km]

% ------------------------------ Altitudes --------------------------------
h_mean_km = a_km - Re;
rp = a_km .* (1 - ecc);
ra = a_km .* (1 + ecc);
hp_km = rp - Re;
ha_km = ra - Re;

% ------------------------------ Output -----------------------------------
tle.norad_id  = norad_id;
tle.epoch     = epoch;
tle.i_deg     = i_deg;
tle.raan_deg  = raan;
tle.e         = ecc;
tle.argp_deg  = argp;
tle.M_deg     = Mdeg;
tle.n_rev_day = nrev;
tle.a_km      = a_km;
tle.h_mean_km = h_mean_km;
tle.hp_km     = hp_km;
tle.ha_km     = ha_km;

end


