%
% Name :
%   ray_test_3d_sp.m
%
% Purpose :
%   Example of using raytrace_3d_sp for a fan of rays. 
%
% Calling sequence :
%   ray_test_3d_sp
%
% Inputs :
%   None
%
% Outputs :
%   None
%
% Modification History:
%   V1.0  M.A. Cervera  22/04/2015
%     Initial version branched from raytrace_3d (V1.2)
%
%   V2.0  M.A. Cervera  03/05/2016
%     Modified to make use of multi-threaded raytrace_3d. IRI2016 is now used
%     to generate the ionosphere.
%



%
% setup general stuff
%
UT = [2000 9 21 0 0];           % UT - year, month, day, hour, minute
speed_of_light = 2.99792458e8;
re = 6376000;                   % radius of Earth to use for coord conv.
R12 = 100;
kp = 3;                         % kp for irregularities model (3 = quiet)
elevs = [3:1:81];               % initial ray elevation
freqs = ones(size(elevs))*15;   % frequency (MHz)
ray_bears = zeros(size(elevs)); % initial bearing of rays
origin_lat = -20.0;             % latitude of the start point of ray
origin_long = 130.0;            % longitude of the start point of ray
origin_ht = 0.0;
doppler_flag = 1;               % interested in Doppler shift

%
% generate ionospheric, geomagnetic and irregularity grids
%
ht_start = 60;          % start height for ionospheric grid (km)
ht_inc = 2;             % height increment (km)
num_ht = 200;           % number of heights (must be <= 401)
lat_start = -20.0;
lat_inc = 0.3;
num_lat = 100.0;        % number of latitudes (must be <= 701)
lon_start= 128.0;
lon_inc = 1.0;
num_lon = 5.0;          % number of longitudes (must be <= 701)
iono_grid_parms = [lat_start, lat_inc, num_lat, lon_start, lon_inc, num_lon, ...
      ht_start, ht_inc, num_ht, ];

B_ht_start = ht_start;          % start height for geomagnetic grid (km)
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

tic
fprintf('Generating ionospheric and geomag grids... ')
[iono_pf_grid, iono_pf_grid_5, collision_freq, Bx, By, Bz] = ...
    gen_iono_grid_3d(UT, R12, iono_grid_parms, ...
                     geomag_grid_parms, doppler_flag, 'iri2016');
toc

% convert plasma frequency grid to  electron density in electrons/cm^3
iono_en_grid = iono_pf_grid.^2 / 80.6164e-6;
iono_en_grid_5 = iono_pf_grid_5.^2 / 80.6164e-6;

%
% call raytrace
%
nhops = 4;                  % number of hops
tol = [1e-7 0.01 25];       % ODE solver tolerance and min max stepsizes
num_elevs = length(elevs);

% Generate the O mode rays
OX_mode = 1;
fprintf('Generating O-mode rays... ')
tic
[ray_data_O, ray_O, ray_state_vec_O] = ...
  raytrace_3d_sp(origin_lat, origin_long, origin_ht, elevs, ray_bears, ...
              freqs, OX_mode, nhops, tol, re, iono_en_grid, iono_en_grid_5, ...
	      collision_freq, iono_grid_parms, Bx, By, Bz, ...
	      geomag_grid_parms);
toc

for rayId=1:num_elevs
  num = length(ray_O(rayId).lat);
  ground_range = zeros(1, num);
  lat = ray_O(rayId).lat;
  lon = ray_O(rayId).lon; 
  ground_range(2:num) = latlon2raz(lat(2:num), lon(2:num), origin_lat, ...
      origin_long, 'spherical', re)/1000.0;
  ray_O(rayId).gndrng = ground_range;
end


% Generate the X mode rays - note in the raytrace_3d call the ionosphere does
% not need to be passed in again as it is already in memory
OX_mode = -1;
 
fprintf('Generating %d X-mode rays ...', num_elevs);
tic
[ray_data_X, ray_X] = ...
  raytrace_3d_sp(origin_lat, origin_long, origin_ht, elevs, ray_bears, ...
                 freqs, OX_mode, nhops, tol, re);
toc
for rayId=1:num_elevs
  num = length(ray_X(rayId).lat);
  ground_range = zeros(1, num);
  lat = ray_X(rayId).lat;
  lon = ray_X(rayId).lon;    
  ground_range(2:num) = latlon2raz(lat(2:num), lon(2:num), origin_lat, ...
      origin_long, 'spherical', re)/1000.0;
  ray_X(rayId).gndrng = ground_range;
end


% Generate the rays for the case where the magnetic field is ignored
OX_mode = 0;

fprintf('Generating ''no-field'' rays... ')

tic
[ray_data_N, ray_N] = ...
  raytrace_3d_sp(origin_lat, origin_long, origin_ht, elevs, ray_bears, ...
                 freqs, OX_mode, nhops, tol, re);
toc
for rayId=1:num_elevs
  num = length(ray_N(rayId).lat);
  ground_range = zeros(1, num);
  lat = ray_N(rayId).lat;
  lon = ray_N(rayId).lon;    
  ground_range(2:num) = latlon2raz(lat(2:num), lon(2:num), origin_lat, ...
      origin_long, 'spherical', re)/1000.0;
  ray_N(rayId).gndrng = ground_range;
end

fprintf('\n')


% finished ray tracing with this ionosphere so clear it out of memory
clear raytrace_3d_sp


% plot the rays
figure(1)
plot3(ray_O(1).lat,  mod(ray_O(1).lon, 360), ray_O(1).height, '.b', ...
      'markersize', 5)
set(gca, 'Zlim', [0 500])
hold on
plot3(ray_X(1).lat, mod(ray_X(1).lon, 360), ray_X(1).height, '.r',  ...
      'markersize',5)
plot3(ray_N(1).lat, mod(ray_N(1).lon, 360), ray_N(1).height, 'g')
for ii = 3:2:num_elevs
  plot3(ray_O(ii).lat, mod(ray_O(ii).lon, 360), ray_O(ii).height,'.b', ...
        'markersize', 5)
  plot3(ray_X(ii).lat, mod(ray_X(ii).lon, 360), ray_X(ii).height, '.r', ...
        'markersize', 5)
  plot3(ray_N(ii).lat, mod(ray_N(ii).lon, 360), ray_N(ii).height, 'g')
end  
hold off
grid on
xlabel('latitude (deg)')
ylabel('longitude (deg)')
zlabel('Height (km)')

figure(2)
start_range = 0;
end_range = 2000;
range_inc = 50;
end_range_idx = fix((end_range-start_range) ./ range_inc) + 1;
start_ht = 0;
start_ht_idx = 1;
height_inc = 5;
end_ht = 350;
end_ht_idx = fix(end_ht ./ height_inc) + 1;
iono_pf_subgrid = zeros(end_ht_idx, end_range_idx);
plot_ray_iono_slice(iono_pf_subgrid, start_range, end_range, range_inc, ...
    start_ht, end_ht, height_inc, ray_O, 'color', [1, 1, 0.99], 'linewidth', 1);
