function [delta_v,delta_v1,delta_v2,VI,VF] = deltaVtot_PLtoAst(t1,t2,id_planet,id_ast)
% This function computes the total delta velocity needed for a transfer
% between a planet and an asteroid.
%
% Inputs: 
%         - t1 [1] [mjd2000]: time of departure
%         - t2 [1] [mjd2000]: time of arrival
%         - id_planet [1]: ephemeris ID of planet of departure
%         - id_ast [1]: ephemeris ID of asteroid of arrival
% Output:
%         - delta_v   [1] [km/s]: total velocity cost of the mission
%         - delta_v1  [1] [km/s]: delta velocity between transfer and departure orbit
%         - delta_v2  [1] [km/s]: delta velocity between arrival and transfer orbit
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

% Calulate orbits through ephemerides
[kep_pl,ksun] = uplanet(t1, id_planet);  % Planet
[kep_a, ~, ~] = ephAsteroids( t2, id_ast ); % Asteroid

% Convert from keplerian to cartesian coordinates
[r_pl, v_pl] = kep2car(kep_pl(1), kep_pl(2),kep_pl(3), kep_pl(4), kep_pl(5), kep_pl(6), ksun);
[r_a, v_a] = kep2car(kep_a(1), kep_a(2),kep_a(3), kep_a(4), kep_a(5), kep_a(6), ksun);

% Lambert problem
orbitType=0; % prograde
Nrev=0;
Ncase=0;
optionsLMR=0;
[~,~,~,~,VI,VF,~,~] = lambertMR(r_pl,r_a,tof,ksun,orbitType,Nrev,Ncase,optionsLMR);

% Compute delta velocity
delta_v1 = norm(VI'-v_pl);
delta_v2 = norm(v_a-VF');
delta_v = delta_v1 + delta_v2;

end
