%
% Name :
%   ray_test_3d_pol.m
%
% Purpose :
%   Example of using raytrace_3d for a single ray. Plots polarization of wave
%   along the ray path
%
% Calling sequence_pol :
%   ray_test_3d
%
% Inputs :
%   None
%
% Outputs :
%   None
%
% Author:
%   V1.0  M.A. Cervera  28/06/2010
%
%   V1.1 M.A. Cervera  19/05/2011
%      More efficient handling of ionospheric and geomagnetic grids grids in
%      call to raytrace_3d  
%
%   V2.0 M.A. Cervera  03/05/2016
%     Modified to make use of multi-threaded raytrace_3d. IRI2016 is now used
%     to generate the ionosphere.
%

%
% setup general stuff
%
UT = [2000 9 21 0 0];        % UT - year, month, day, hour, minute
speed_of_light = 2.99792458e8;
R12 = 100;
kp = 3;                      % kp for irregularities model (3 = quiet)
elev = 25;              % initial ray elevation - two rays
freq = 10.0;          % frequency (MHz)
ray_bear = 0;            % initial bearing of radar ray


origin_lat = -20.0;          % latitude of the start point of ray
origin_long = 130.0;         % longitude of the start point of ray
origin_ht = 0.0;

doppler_flag = 0;            % not interested in Doppler shift and spread

%
% generate ionospheric, geomagnetic and irregularity grids
%
ht_start = 60;          % start height for ionospheric grid (km)
ht_inc = 2;             % height increment (km)
num_ht = 200;           % number of  heights (must be <= 401)

lat_start = -22.0;
lat_inc = 0.5;
num_lat = 100.0;        % number of latitudes (must be <= 701)

lon_start= 128.0;
lon_inc = 1.0;
num_lon = 5.0;          % number of longitudes (must be <= 701)
 
iono_grid_parms = [lat_start, lat_inc, num_lat, lon_start, lon_inc, num_lon, ...
      ht_start, ht_inc, num_ht, ];

B_ht_start = 0; %ht_start;          % start height for geomagnetic grid (km)
B_ht_inc = 10;                  % height increment (km)
B_num_ht = ceil(num_ht .* ht_inc ./ B_ht_inc) + 1;
B_lat_start = lat_start;
B_lat_inc = 1.0;
B_num_lat = ceil(num_lat .* lat_inc ./ B_lat_inc) + 1;
B_lon_start = lon_start;
B_lon_inc = 1.0;
B_num_lon = ceil(num_lon .* lon_inc ./ B_lon_inc) + 1; 
geomag_grid_parms = [B_lat_start, B_lat_inc, B_num_lat, B_lon_start, ...
      B_lon_inc, B_num_lon, B_ht_start, B_ht_inc, B_num_ht];

indices = [R12, 0, 0, 0, kp];
doppler_flag = 0;

tic
fprintf('Generating ionospheric and geomag grids... ')
[iono_pf_grid, iono_pf_grid_5, collision_freq, Bx, By, Bz] = ...
    gen_iono_grid_3d(UT, R12, iono_grid_parms, geomag_grid_parms, ...
                     doppler_flag, 'iri2016');
toc

% convert plasma frequency grid to  electron density in electrons/cm^3
iono_en_grid = iono_pf_grid.^2 / 80.6164e-6;
iono_en_grid_5 = iono_pf_grid_5.^2 / 80.6164e-6;



%
% call raytrace
%
nhops = 1;           % number of hops
tol = [1e-7 0.01 10];          % rkf tolerance

ray_O = [];
ray_X = [];
ray_N = [];

%
% Generate O mode ray
%
OX_mode = 1;
fprintf('Generating O-mode rays... ')
tic

% first call to raytrace_3d so pass in the ionosphere
[ray_O, ray_path_O, ray_state_vec_O] = ...
    raytrace_3d(origin_lat, origin_long, origin_ht, elev, ray_bear, freq, ...
	     OX_mode, nhops, tol, iono_en_grid, iono_en_grid_5, ...
	     collision_freq, iono_grid_parms, Bx, By, Bz, geomag_grid_parms); 
toc


%
% Generate X mode ray
%
OX_mode = -1;
fprintf('Generating X-mode rays... ')

% ionosphere already in memory so no need to pass it in again
[ray_X, ray_path_X, ray_state_vec_X] = ...
    raytrace_3d(origin_lat, origin_long, origin_ht, elev, ray_bear, freq, ...
	     OX_mode, nhops, tol);
toc


%
% Generate 'no-field' mode ray
%
OX_mode = 0;
fprintf('Generating ''no-field'' rays... ')
tic

% ionosphere already in memory so no need to pass it in again
[ray_N, ray_path_N, ray_state_vec_N] = ...
    raytrace_3d(origin_lat, origin_long, origin_ht, elev, ray_bear, freq, ...
	     OX_mode, nhops, tol);

toc
fprintf('\n')

% ionosphere is no longer needed by ratrace_3d so clear it.
clear raytrace_3d


%
% plot the rays
%
figure(1)
plot3(ray_path_O(1).lat, ray_path_O(1).lon, ray_path_O(1).height, '.', ...
      'markersize', 5)
hold on
plot3(ray_path_X(1).lat, ray_path_X(1).lon, ray_path_X(1).height, '.r',  ...
      'markersize',5)
plot3(ray_path_N(1).lat, ray_path_N(1).lon, ray_path_N(1).height, '.g')
hold off
grid on
xlabel('latitude (deg)')
ylabel('longitude (deg)')
zlabel('Height (km)')
legend('O-mode', 'X-mode', 'No-field')


figure(2)
set(gcf, 'Position', [20 420 850 650])
subplot(2,2,1)
plot(ray_path_O(1).group_range, abs(ray_path_O(1).polariz_mag), '.')
xlabel('Group Range (km)')
ylabel('Polarization, |R|')
title('Wave {\bf E}-field vector axial ratio')
grid on
hold on 
plot(ray_path_X(1).group_range, abs(ray_path_X(1).polariz_mag), 'r.')
hold off
legend('O-mode', 'X-mode')

subplot(2,2,2)
plot(ray_path_O(1).group_range, ray_path_O(1).wave_Efield_tilt, '.')
xlabel('Group Range (km)')
ylabel('Polarization, \psi (degrees)')
title('Wave {\bf E}-field forward tilt angle') 
grid on
hold on 
plot(ray_path_X(1).group_range, ray_path_X(1).wave_Efield_tilt, 'r.')
hold off
legend('O-mode', 'X-mode')

subplot(2,2,3)
plot(ray_path_O(1).group_range, ray_path_O(1).wavenorm_B_angle, '.')
title('Angle between wave-normal and {\bf B}-field') 
xlabel('Group Range (km)')
ylabel('\theta (degrees)')
grid on
hold on 
plot(ray_path_X(1).group_range, ray_path_X(1).wavenorm_B_angle, 'r.')
hold off
legend('O-mode', 'X-mode')

subplot(2,2,4)
plot(ray_path_O(1).group_range, ray_path_O(1).wavenorm_ray_angle, '.')
title('Angle between wave-normal and ray direction') 
xlabel('Group Range (km)')
ylabel('\alpha (degrees)')
grid on
hold on 
plot(ray_path_X(1).group_range, ray_path_X(1).wavenorm_ray_angle, 'r.')
hold off
legend('O-mode', 'X-mode')


figure(3)
set(gcf, 'Position', [400 600 850 325])
subplot(1,2,1)
plot(ray_path_O(1).group_range, ray_path_O(1).refractive_index, '.')
xlabel('Group Range (km)')
ylabel('Refractive Index')
title('Refractive Index')
grid on
hold on 
plot(ray_path_X(1).group_range, ray_path_X(1).refractive_index, 'r.')
hold off
legend('O-mode', 'X-mode')

subplot(1,2,2)
plot(ray_path_O(1).group_range, ray_path_O(1).group_refractive_index, '.')
xlabel('Group Range (km)')
ylabel('Group Refractive Index')
title('Group Refractive Index') 
grid on
hold on 
plot(ray_path_X(1).group_range, ray_path_X(1).group_refractive_index, 'r.')
hold off
legend('O-mode', 'X-mode')
