function [delta_v,delta_v1,delta_v2,VI,VF] = deltaVtot_pl2pl(t1,t2,planet_dep, planet_arr)
% This function computes the total delta velocity needed for a transfer
% between two planets of the Solar system.
%
% Inputs: 
%         - t1 [1] [mjd2000]: time of departure
%         - t2 [1] [mjd2000]: time of arrival
%         - planet_dep [1]: ephemeris ID of planet of departure 
%         - planet_arr [1]: ephemeris ID of planet of arrival 
% Output:
%         - delta_v  [1] [km/s]: total velocity cost of the mission
%         - delta_v1 [1] [km/s]: delta velocity between transfer and departure orbit
%         - delta_v2 [1] [km/s]: delta velocity between arrival and transfer orbit
%         - VI [3] [km/s]: initial velocity of transfer arc in cartesian coordinates
%         - VF [3] [km/s]: final velocity of transfer arc in cartesian coordinates
%
% CONTRIBUTORS:
% -Simone Frescura
% -Lorenzo Veronese
% -Nicola Pifano
% -Luca Puzzi
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Compute time of flight 
tof = (t2 - t1) * 86400; % [s]

% Calculate orbits through ephemerides
[kep1,~] = uplanet(t1, planet_dep);  
[kep2,ksun] = uplanet(t2, planet_arr); 

% Convert from keplerian to cartesian coordinates
[r1, v1] = kep2car(kep1(1), kep1(2),kep1(3), kep1(4), kep1(5), kep1(6), ksun);
[r2, v2] = kep2car(kep2(1), kep2(2), kep2(3), kep2(4), kep2(5), kep2(6), ksun);

% Lambert problem
orbitType=0; % prograde
Nrev=0;
Ncase=0;
optionsLMR=0;
[~,~,~,~,VI,VF,~,~] = lambertMR(r1,r2,tof,ksun,orbitType,Nrev,Ncase,optionsLMR);

% Compute delta velocity
delta_v1 = norm(VI'-v1);
delta_v2 = norm(v2-VF');
delta_v = delta_v1 + delta_v2;

end
