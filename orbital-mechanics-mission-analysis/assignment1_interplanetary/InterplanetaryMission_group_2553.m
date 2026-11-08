
% Main script for Assignment 1 – Interplanetary Explorer Mission
%
% Group 2553
%
% AUTHORS:
% - Simone Frescura
% - Lorenzo Veronese 
% - Nicola Pifano
% - Luca Puzzi
%

clear;     
close all;
clc;        

% Project functions and data (run this script from its own folder);
% the course-provided functions go in ../course_functions (see README)
addpath('functions', fullfile('..','course_functions'));

%% 1. Celestial objects
id_ven  = 2;        % Departure planet: VENUS
id_mar  = 4;        % Fly-by planet: MARS
id_ast  = 184430;   % Asteroid's ID

%% 2. Costants

muSun   = astroConstants(4);   % [km^3/s^2] 
muVenus = astroConstants(12);   % [km^3/s^2] 
muMars  = astroConstants(14);   % [km^3/s^2] 

rVenus  = astroConstants(22);   % [km] radius of Venus
rMars   = astroConstants(24);   % [km] radius of Mars

hmin = 230;  %[km] minimum altitude of fly-by

%% 3. Time window
% Format: [year month day hour minuts seconds]
dep_ven_i = [2030 1 1 0 0 0];
dep_ven_f = [2060 1 1 0 0 0];
fb_mar_i = [2030 1 2 0 0 0];
fb_mar_f = [2060 1 1 0 0 0];
arr_ast_i = [2030 1 3 0 0 0];
arr_ast_f = [2060 1 1 0 0 0];

t_ven_min=date2mjd2000(dep_ven_i);
t_ven_max=date2mjd2000(dep_ven_f);
t_mar_min=date2mjd2000(fb_mar_i);
t_mar_max=date2mjd2000(fb_mar_f);
t_ast_min=date2mjd2000(arr_ast_i);
t_ast_max=date2mjd2000(arr_ast_f);

%% Coarse grid search to define better time window

% Grid resolution
N_dep=20; 
N_fb=20;
N_arr=20;

tdep_vec=linspace(t_ven_min,t_ven_max,N_dep);
tfb_vec=linspace(t_mar_min,t_mar_max,N_fb);
tarr_vec=linspace(t_ast_min,t_ast_max,N_arr);

DV_map = NaN(N_dep,N_fb,N_arr); 

for i=1:N_dep
  for j=1:N_fb
    for k=1:N_arr
      t = [tdep_vec(i); tfb_vec(j); tarr_vec(k)];
       if (t(2) <= t(1)) || (t(3) <= t(2))
                DV_map(i,j,k) = Inf;
                continue
       end
      DV_map(i,j,k) = mission_objective(t,id_ven,id_mar,id_ast,muMars,rMars,hmin);
    end
  end
end

DV_lim = 25; % choose a limit deltav
good = DV_map<DV_lim; % select deltav respecting the constraint

% Departure window
idx_dep = any(any(good,2),3);
tdep_good = tdep_vec(idx_dep);
tdep_min=min(tdep_good); % lower boundary of departure window
tdep_max=max(tdep_good); % upper boundary of departure window
fprintf('Departure window:\n %d %d %d %d %d %f\n%d %d %d %d %d %f\n', ...
    mjd20002date(tdep_min),mjd20002date(tdep_max));

% Fly-by window
idx_fb = any(any(good,1),3);
tfb_good = tfb_vec(idx_fb);
tfb_min=min(tfb_good); % lower boundary of fly-by window
tfb_max=max(tfb_good); % upper boundary of fly-by window
fprintf('Fly-by window:\n %d %d %d %d %d %f\n%d %d %d %d %d %f\n', ...
    mjd20002date(tfb_min),mjd20002date(tfb_max));

% Arrival window
idx_arr = any(any(good,1),2);
tarr_good = tarr_vec(idx_arr);
tarr_min=min(tarr_good); % lower boundary of arrival window
tarr_max=max(tarr_good); % upper boundary of arrival window
fprintf('Arrival window:\n %d %d %d %d %d %f\n%d %d %d %d %d %f\n', ...
    mjd20002date(tarr_min),mjd20002date(tarr_max));

lb = [tdep_min; tfb_min; tarr_min]; % lower boundary
ub = [tdep_max; tfb_max; tarr_max]; % upper boundary

%% Grid search to find best initial condition
DVbest = Inf;
N_dep=30; 
N_fb=30;
N_arr=30;
tspan_dep=linspace(min(tdep_good),max(tdep_good),N_dep);
tspan_fb=linspace(min(tfb_good),max(tfb_good),N_fb);
tspan_arr=linspace(min(tarr_good),max(tarr_good),N_arr);

for t1 = tspan_dep
    for t2 = tspan_fb
        for t3 = tspan_arr
            x0 = [t1; t2; t3];
            DV = mission_objective(x0,id_ven,id_mar,id_ast,muMars,rMars,hmin);
            if DV < DVbest
                DVbest = DV;
                x0_best = x0;
            end
        end
    end
end

fprintf('Date of departure from Venus: %d %d %d %d %d %f\n',mjd20002date(x0_best(1)))
fprintf('Date of fly-by: %d %d %d %d %d %f\n',mjd20002date(x0_best(2)))
fprintf('Date of arrival on the asteroid: %d %d %d %d %d %f\n',mjd20002date(x0_best(3)))
fprintf('Delta v of the mission: %f km/s\n',DVbest)

% Trial of different initial conditions
%x0_best=[(ub(1)+lb(1))/2,(ub(2)+lb(2))/2,(ub(3)+lb(3))/2]; % --> deltaV_tot=15.5 km/s
%x0_best=[ub(1),ub(2),ub(3)]; % --> deltaV_tot=17.8 km/s
%x0_best=[lb(1)+300,lb(2)+300,lb(3)+300]; % --> deltaV_tot=31.4 km/s


%% 4. Function to minimize


options = optimoptions('fmincon','Algorithm','sqp','Display','iter','MaxFunctionEvaluations',3e4);


[x_opt, deltav_min] = fmincon(@(x)mission_objective(x,id_ven,id_mar,id_ast,muMars,rMars,hmin),x0_best, ...
                        [], [], [], [], lb, ub,[],options);

%% 5. Results
t_dep = x_opt(1);
t_fb  = x_opt(2);
t_arr = x_opt(3);
tof_mission=(t_arr-t_dep)/365; %[years]

fprintf('Date of departure from Venus: %d %d %d %d %d %f\n',mjd20002date(t_dep))
fprintf('Date of fly-by: %d %d %d %d %d %f\n',mjd20002date(t_fb))
fprintf('Date of arrival on the asteroid: %d %d %d %d %d %f\n',mjd20002date(t_arr))
fprintf('Delta v of the mission: %f km/s\n',deltav_min)
fprintf('Time of flight of the mission: %f years\n',tof_mission)

%%  Venus orbit
% Venus at departure
[kep_ven_dep,~] = uplanet(t_dep, id_ven);
[r_ven_dep,v_ven_dep] = kep2car(kep_ven_dep(1),kep_ven_dep(2), ...
         kep_ven_dep(3),kep_ven_dep(4),kep_ven_dep(5),kep_ven_dep(6),muSun);
% Venus at fly-by
[kep_ven_fb,~] = uplanet(t_fb, id_ven);
[r_ven_fb,~] = kep2car(kep_ven_fb(1),kep_ven_fb(2), ...
         kep_ven_fb(3),kep_ven_fb(4),kep_ven_fb(5),kep_ven_fb(6),muSun);
% Venus at arrival
[kep_ven_arr,~] = uplanet(t_arr, id_ven);
[r_ven_arr,~] = kep2car(kep_ven_arr(1),kep_ven_arr(2), ...
         kep_ven_arr(3),kep_ven_arr(4),kep_ven_arr(5),kep_ven_arr(6),muSun);

y0_ven=[r_ven_dep; v_ven_dep];
T_ven=2*pi*sqrt( kep_ven_dep(1)^3/muSun ); % orbital period
options = odeset('RelTol',1e-13,'AbsTol',1e-14);
tspan_ven=linspace(0,T_ven,1000);
[time_ven,Y_ven]=ode113(@(t,y) ode_2bp(t,y,muSun),tspan_ven,y0_ven,options);

%% Mars orbit
% Mars at departure
[kep_mar_dep,~] = uplanet(t_dep, id_mar);
[r_mar_dep,~] = kep2car(kep_mar_dep(1),kep_mar_dep(2),kep_mar_dep(3), ...
    kep_mar_dep(4),kep_mar_dep(5),kep_mar_dep(6),muSun);

% Mars at flyby
[kep_mar_fb,~] = uplanet(t_fb, id_mar);
[r_mar_fb,v_mar_fb] = kep2car(kep_mar_fb(1),kep_mar_fb(2),kep_mar_fb(3), ...
    kep_mar_fb(4),kep_mar_fb(5),kep_mar_fb(6),muSun);

% Mars at arrival
[kep_mar_arr,~] = uplanet(t_arr, id_mar);
[r_mar_arr,~] = kep2car(kep_mar_arr(1),kep_mar_arr(2),kep_mar_arr(3), ...
    kep_mar_arr(4),kep_mar_arr(5),kep_mar_arr(6),muSun);


y0_mar=[r_mar_fb; v_mar_fb];
T_mar=2*pi*sqrt( kep_mar_dep(1)^3/muSun );

tspan_mar=linspace(0,T_mar,1000);
[time_mar,Y_mar]=ode113(@(t,y) ode_2bp(t,y,muSun),tspan_mar,y0_mar,options);

%% Asteroid orbit
% Asteroid at departure
[kep_ast_dep,~] = ephAsteroids(t_dep, id_ast);
[r_ast_dep,~] = kep2car(kep_ast_dep(1),kep_ast_dep(2),kep_ast_dep(3), ...
    kep_ast_dep(4),kep_ast_dep(5),kep_ast_dep(6),muSun);

% Asteroid at flyby
[kep_ast_fb,~] = ephAsteroids(t_fb, id_ast);
[r_ast_fb,~] = kep2car(kep_ast_fb(1),kep_ast_fb(2),kep_ast_fb(3), ...
    kep_ast_fb(4),kep_ast_fb(5),kep_ast_fb(6),muSun);

% Asteroid at arrival
[kep_ast_arr,~] = ephAsteroids(t_arr, id_ast);
[r_ast_arr,v_ast_arr] = kep2car(kep_ast_arr(1),kep_ast_arr(2),kep_ast_arr(3), ...
    kep_ast_arr(4),kep_ast_arr(5),kep_ast_arr(6),muSun);

y0_ast=[r_ast_arr; v_ast_arr];
T_ast=2*pi*sqrt( kep_ast_arr(1)^3/muSun );

tspan_ast=linspace(0,T_ast,1000);
[time_ast,Y_ast]=ode113(@(t,y) ode_2bp(t,y,muSun),tspan_ast,y0_ast,options);

%% Plot of starting orbits
figure(1)
plot3(Y_ven(:,1),Y_ven(:,2),Y_ven(:,3)) % Venus orbit
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Planets orbit');
axis equal;
grid on;
hold on;
plot3(Y_mar(:,1),Y_mar(:,2),Y_mar(:,3)) % Mars orbit
plot3(Y_ast(:,1),Y_ast(:,2),Y_ast(:,3)) % Asteroid orbit
plot3(0,0,0,'o','MarkerSize', 10,'MarkerEdgeColor', 'y', 'MarkerFaceColor', 'y') % Sun
legend('Venus','Mars','Asteroid','Sun')

%% First lambert arc
orbitType=0; % prograde
Nrev=0;
Ncase=0;
optionsLMR=0;
tof1 = (t_fb - t_dep) * 86400;
[~,~,~,~,VI1,VF1,~,transfer_angle1] = lambertMR(r_ven_dep,r_mar_fb,tof1,muSun,orbitType,Nrev,Ncase,optionsLMR);
y0_tr1=[r_ven_dep; VI1'];

[kep_arc1(1),kep_arc1(2),kep_arc1(3),kep_arc1(4),kep_arc1(5),~]=car2kep(r_ven_dep,VI1',muSun);
kep_arc1(3)=rad2deg(kep_arc1(3));
kep_arc1(4)=rad2deg(kep_arc1(4));
kep_arc1(5)=rad2deg(kep_arc1(5));

tspan_tr1=linspace(0,tof1,1000);
[time_tr1,Y_tr1]=ode113(@(t,y) ode_2bp(t,y,muSun),tspan_tr1,y0_tr1,options);


%% Second lambert arc
tof2 = (t_arr - t_fb) * 86400;
[~,~,~,~,VI2,VF2,~,transfer_angle2] = lambertMR(r_mar_fb,r_ast_arr,tof2,muSun,orbitType,Nrev,Ncase,optionsLMR);
y0_tr2=[r_mar_fb; VI2'];

[kep_arc2(1),kep_arc2(2),kep_arc2(3),kep_arc2(4),kep_arc2(5),kep_arc2(6)]=car2kep(r_mar_fb,VI2',muSun);
kep_arc2(3)=rad2deg(kep_arc2(3));
kep_arc2(4)=rad2deg(kep_arc2(4));
kep_arc2(5)=rad2deg(kep_arc2(5));

tspan_tr2=linspace(0,tof2,1000);
[time_tr2,Y_tr2]=ode113(@(t,y) ode_2bp(t,y,muSun),tspan_tr2,y0_tr2,options);

%% Plot Heliocentric arc
figure(2)
plot3(Y_ven(:,1),Y_ven(:,2),Y_ven(:,3)) % Venus orbit
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Heliocentric orbits');
axis equal;
grid on;
hold on;
plot3(Y_mar(:,1),Y_mar(:,2),Y_mar(:,3)) % Mars orbit
plot3(Y_ast(:,1),Y_ast(:,2),Y_ast(:,3)) % Asteroid orbit
plot3(0,0,0,'o','MarkerSize', 10,'MarkerEdgeColor', 'y', 'MarkerFaceColor', 'y') % Sun
plot3(Y_tr1(:,1),Y_tr1(:,2),Y_tr1(:,3)) % First transfer orbit
plot3(Y_tr2(:,1),Y_tr2(:,2),Y_tr2(:,3)) % Second transfer orbit

plot3(r_ven_dep(1), r_ven_dep(2), r_ven_dep(3), ...
    'o','MarkerSize', 8, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'k') % venus at departure
plot3(r_ven_fb(1), r_ven_fb(2), r_ven_fb(3), ...
    'o','MarkerSize', 8, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'k') % venus at fly-by
plot3(r_ven_arr(1), r_ven_arr(2), r_ven_arr(3), ...
    'o','MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k') % venus at arrival

plot3(r_mar_dep(1), r_mar_dep(2), r_mar_dep(3), ...
    'o','MarkerSize', 8, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'k') % Mars at dep
plot3(r_mar_fb(1), r_mar_fb(2), r_mar_fb(3), ...
    'o','MarkerSize', 8, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'k') % Mars at fly-by
plot3(r_mar_arr(1), r_mar_arr(2), r_mar_arr(3), ...
    'o','MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k') % Mars at arrival

plot3(r_ast_dep(1), r_ast_dep(2), r_ast_dep(3), ...
    'o', 'MarkerSize', 8, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'k') % Asteroid at departure
plot3(r_ast_fb(1), r_ast_fb(2), r_ast_fb(3), ...
    'o', 'MarkerSize', 8, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'k') % Asteroid at fly-by
plot3(r_ast_arr(1), r_ast_arr(2), r_ast_arr(3), ...
    'o', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k') % Asteroid at arrival

legend('Venus','Mars','Asteroid','Sun','First transfer orbit',' Second transfer orbit', ...
    'planet at departure time','planet at fly-by time','planet at arrival time')


%% FLY-BY

orbit_rad_M= 228*10^6;       % [km] mean orbital radius of Mars
r_SOI = orbit_rad_M*0.0025;  % [km] radius of SOI (0.0025 taken from lab4)

V_minus=VF1';  % Incoming heliocentric velocity
V_plus=VI2';   % Outcoming heliocentric velocity

v_inf_minus = V_minus - v_mar_fb;  % Velocity at infinite before fly-by
v_inf_plus = V_plus - v_mar_fb;    % Velocity at infinite after fly-by

%  Turn angle 
delta = acos(dot(v_inf_minus,v_inf_plus)/(norm(v_inf_plus)*norm(v_inf_minus)));

% Radius of pericentre
rp_min = rMars + hmin;
fun = @(rp) (2*asin(1/(1 + (rp*norm(v_inf_minus)^2)/muMars)))/2 + ...
      (2*asin(1/(1 + (rp*norm(v_inf_plus)^2)/muMars)))/2 - delta;  % delta_m/2 + delta_p/2 - delta =0
r_p=fzero(fun,rp_min + 1000);
fprintf('Radius of pericentre: %f km\n',r_p);

h_fb=r_p-rMars;  % [km] Altitude of fly-by
fprintf('Altitude of fly-by: %f km\n',h_fb)

% delta velocity at pericentre
vp_minus=sqrt(norm(v_inf_minus)^2 + 2*muMars/r_p); % [km/s] velocity at pericentre (minus)
vp_plus=sqrt(norm(v_inf_plus)^2 + 2*muMars/r_p); % [km/s] velocity at pericentre (plus)
delta_vp=abs(vp_plus - vp_minus); % [km/s] delta velocity of maneuvre at pericentre 

deltaV_fb=norm(V_plus-V_minus); % delta velocity due to fly-by

fprintf('Velocity change due to the fly-by: %f km/s\n',deltaV_fb)
fprintf('Velocity change due to powered maneuvre at pericentre: %e km/s\n',delta_vp)

% Incoming hyperbola
a_minus = -muMars/(norm(v_inf_minus)^2);       % semi-major axis 
e_minus = 1 + (r_p*norm(v_inf_minus)^2)/muMars;  % eccentricity
theta_inf_minus = acos(-1/e_minus);
theta_minus = linspace(-theta_inf_minus,0, 1000); % [rad] Angular range
r_hyp_minus = a_minus*(1 - e_minus^2) ./ (1 + e_minus*cos(theta_minus)); % [km] hyperbolic radius
x_hyp_minus = r_hyp_minus .* cos(theta_minus); % Hyperbolic x coordinate
y_hyp_minus = r_hyp_minus .* sin(theta_minus); % Hyperbolic y coordinate

% Outocoming hyperbola
a_plus = -muMars/(norm(v_inf_plus)^2);         % semi-major axis  
e_plus = 1 + (r_p*norm(v_inf_plus)^2)/muMars;  % eccentricity
theta_inf_plus = acos(-1/e_plus);
theta_plus = linspace(0, theta_inf_plus, 1000); % [rad] Angular range
r_hyp_plus = a_plus*(1 - e_plus^2) ./ (1 + e_plus*cos(theta_plus)); % [km] hyperbolic radius
x_hyp_plus = r_hyp_plus .* cos(theta_plus);  % Hyperbolic x coordinate
y_hyp_plus = r_hyp_plus .* sin(theta_plus);  % Hyperbolic y coordinate

% Duration of fly-by
E_minus = acosh( (r_SOI/abs(a_minus) + 1)/e_minus );  % [rad] Hyperbolic anomaly of entrance point in SOI
t_minus = abs(sqrt((abs(a_minus))^3/muMars) * ( e_minus*sinh(E_minus) - E_minus )); % [s]

E_plus = acosh( (r_SOI/abs(a_minus) + 1)/e_plus );  % [rad] Hyperbolic anomaly of exit point of SOI
t_plus = abs(sqrt((abs(a_plus))^3/muMars) * ( e_plus*sinh(E_plus) - E_plus )); % [s]

tof_fb=t_plus+t_minus;  % [s] time duration of fly-by
fprintf('Time duration of the fly-by: %f hours\n',tof_fb/3600);

% Delta velocity at departure
deltaV_dep= norm(VI1'-v_ven_dep);
fprintf('Delta velocity at departure: %f km/s\n',deltaV_dep);

% Delta velocity at arrival
deltaV_arr=norm(v_ast_arr-VF2');
fprintf('Delta velocity at arrival: %f km/s\n',deltaV_arr);


figure(3)
marsImg = imread(fullfile('data','2k_mars.jpg'));
% Sphere mesh
[Xe,Ye,Ze] = sphere(60);
Xe_rot = Xe;
Ye_rot = Ye;
Ze_rot = -Ze;
h = surf(rMars*Xe_rot, rMars*Ye_rot, rMars*Ze_rot);
set(h, 'FaceColor', 'texturemap', 'CData', marsImg,'EdgeColor', 'none');


hold on;
axis equal;
grid on;
plot(x_hyp_minus(2:end), y_hyp_minus(2:end), 'LineWidth', 2); 
hold on;
plot(x_hyp_minus(end), y_hyp_minus(end), 'o');
plot(x_hyp_plus(1:end-1), y_hyp_plus(1:end-1),'LineWidth', 2);
xlabel('x [km]');
ylabel('y [km]');
zlabel('z [km]');
title('Hyperbolic arcs');
legend('Mars','Incoming hyperbola','Pericentre','Outocoming hyperbola')
xlim([-30000 30000])
ylim([-30000 30000])

%% Numerical results to send
fprintf('NUMERICAL RESULTS TO SEND:\n')
fprintf('Date of departure in MJD2000: %f\n',t_dep)
fprintf('Date of fly-by in MJD2000: %f\n',t_fb)
fprintf('Date of arrival in MJD2000: %f\n',t_arr)
fprintf('Delta v of the mission: %f km/s\n',deltav_min)
fprintf('Radius of pericentre: %f km\n',r_p);
fprintf('Time duration of the fly-by: %f hours\n',tof_fb/3600);
ratio = delta_vp/deltav_min;
fprintf('Ratio between deltav of impulsive manoeuvre and total deltav: %f',ratio)