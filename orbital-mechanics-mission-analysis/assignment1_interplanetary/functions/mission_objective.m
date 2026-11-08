function delta_vtot = mission_objective(t,id_ven,id_mar,id_ast,mu_M,R_M,h_min)
% This function gives the total delta velocity needed to transfer from Venus 
% to asteroid, by doing a powered fly-by on Mars.
%
% INPUTS: 
%        - t [3][mjd2000]: vector of time of departure, time of fly-by and
%                          time of arrival
%        - id_ven [1]: ephemerides id of Venus
%        - id_mar [1]: ephemerides id of Mars
%        - id_ast [1]: ephemerides id of asteroid
%        - mu_M [1][km^3/s^2]: planetary constant of Mars
%        - R_M [1][km]: mean radius of Mars
%        - h_min [1]: minimum altitude of fly-by
%
% OUTPUT:
%        - delta_vtot [1][km/s]: total delta velocity
%
% CONTRIBUTORS:
% -Simone Frescura
% -Lorenzo Veronese
% -Nicola Pifano
% -Luca Puzzi
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

t1 = t(1);
t2 = t(2);
t3 = t(3);

% First heliocentric arc
[~,delta_v_arc1,~,~,VF_m] =deltaVtot_pl2pl(t1,t2,id_ven,id_mar);

% Second heliocentri arc
[~,~,delta_v_arc2,VI_m,~] = deltaVtot_PLtoAst(t2,t3,id_mar,id_ast);

% Orbit and speed of Mars
[kep_mars,mu_Sun] = uplanet(t2, id_mar);  
[~, V_mars] = kep2car(kep_mars(1), kep_mars(2),kep_mars(3), kep_mars(4), kep_mars(5), kep_mars(6), mu_Sun);

% hyperbolic speeds
v_inf_minus = VF_m - V_mars';
v_inf_plus  = VI_m - V_mars';

% Turn angle 
delta = acos(dot(v_inf_minus,v_inf_plus)/(norm(v_inf_plus)*norm(v_inf_minus)));

% Feasibility check
rp_min = R_M + h_min;

delta_max = asin(1/(1 + rp_min*norm(v_inf_minus)^2/mu_M)) + asin(1/(1 + rp_min*norm(v_inf_plus)^2/mu_M));

% IMPOSSIBLE FLYBY 
if delta > delta_max || ~isreal(delta)
    delta_vtot = 1e5 + 1e3*(delta - delta_max)^2;
    return
end

% Radius of pericentre
fun = @(rp) (2*asin(1/(1 + (rp*norm(v_inf_minus)^2)/mu_M)))/2 + ...
      (2*asin(1/(1 + (rp*norm(v_inf_plus)^2)/mu_M)))/2 - delta;  % delta_m/2 + delta_p/2 - delta =0
r_p=fzero(fun,rp_min + 1000);

% Velocities at pericentre
vp_minus = sqrt(norm(v_inf_minus)^2 + 2*mu_M/r_p);
vp_plus  = sqrt(norm(v_inf_plus)^2 + 2*mu_M/r_p);

% Delta v of powered maneuvre
delta_vp = abs(vp_plus-vp_minus);

% Total delta v
delta_vtot = delta_v_arc1 + delta_v_arc2 + delta_vp;


end


