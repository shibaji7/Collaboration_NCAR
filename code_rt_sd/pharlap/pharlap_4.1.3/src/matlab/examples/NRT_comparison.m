%
% Name :
%   NRT_comparison.m
%
% Purpose :
%   Compares the WGS84 2D NRT and 3D (no magnetic field) NRT for a fan of rays
%
% Calling sequence :
%   NRT_comparison
%
% Inputs :
%   None
%
% Outputs :
%   None
%
% Change log:
%   V1.0  M.A. Cervera  06/07/2016
%     Initial Version
%


%
% setup general stuff
%
UT = [2000 9 21 0 0];           % UT - year, month, day, hour, minute
speed_of_light = 2.99792458e8;
R12 = 100;
kp = 3;                              % kp for irregularities model (3 = quiet)
elevs = [3:0.1:81];                  % initial elevation of rays
freqs = ones(size(elevs))*15;        % frequency (MHz)
ray_bears = 30 + zeros(size(elevs)); % initial bearing of rays

origin_lat = -20.0;                  % latitude of the start point of rays
origin_long = 130.0;                 % longitude of the start point of rays
origin_ht = 0.0;                     % altitude of the start point of rays
doppler_flag = 1;                    % interested in Doppler shift


%
% generate 3D ionospheric, geomagnetic and irregularity grids
%
ht_start = 60;          % start height for ionospheric grid (km)
ht_inc = 1;             % height increment (km)
num_ht = 300;           % number of  heights (must be < 201)
lat_start = -21.0;
lat_inc = 0.2;
num_lat = 100.0;
lon_start= 129.0;
lon_inc = 0.2;
num_lon = 100;

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

[iono, iono_extra] = iri2016(lat_start, lon_start, R12, UT, ht_start, ...
    ht_inc,  num_ht);
en_profile = iono(1, :) / 1e6;      % convert to electrons per cm^3

UT5 = UT; UT5(5) = UT5(5) + 5;      % UT 5 minutes later
[iono_5, iono_extra_5] = iri2016(lat_start, lon_start, R12, UT5, ...
                                       ht_start, ht_inc, num_ht);
en_profile_5 = iono_5(1, :) / 1e6;  % convert to electrons per cm^3

% simple empirical formula for collision frequency
NmE = iono_extra(5);
foE = sqrt(NmE  * 80.6164e-12);
height_arr = [ht_start : ht_inc : ht_start + (num_ht-1)*ht_inc];
cf_profile = 4 .* foE.^4 + 1000 .* exp(-(height_arr - 140) ./ 20);

iono_en_grid = zeros(num_lat, num_lon, num_ht);
iono_en_grid_5 = zeros(num_lat, num_lon, num_ht);
collision_freq  = zeros(num_lat, num_lon, num_ht);
for ii = 1:num_lat
  for jj = 1:num_lon
    iono_en_grid(ii, jj, :) = en_profile;
    iono_en_grid_5(ii, jj, :) = en_profile_5;
    collision_freq(ii, jj, :) = cf_profile;
  end
end

% set B field to zero as we are only doing 'no-field' 3D NRT
Bx = zeros(B_num_lat, B_num_lon, B_num_ht);
By = zeros(B_num_lat, B_num_lon, B_num_ht);
Bz = zeros(B_num_lat, B_num_lon, B_num_ht);


%
% call no-field 3D raytrace
%
nhops = 1;                  % number of hops
tol = [1e-9 0.001 25];       % ODE solver tolerance and min max stepsizes
num_elevs = length(elevs);
OX_mode = 0;

fprintf('3D NRT: generating %d ''no-field'' rays ...', num_elevs);
tic
[ray_data_N, ray_N] = ...
  raytrace_3d(origin_lat, origin_long, origin_ht, elevs, ray_bears, ...
              freqs, OX_mode, nhops, tol, iono_en_grid, ...
	      iono_en_grid_5, collision_freq, iono_grid_parms, Bx, By, Bz, ...
	      geomag_grid_parms);
toc
fprintf('\n')

% finished ray tracing with this ionosphere so clear it out of memory
clear raytrace_3d

% plot ground range vs elevation
figure(1)
idx_goodray = find([ray_data_N(:).ray_label] == 1);
group_range = [ray_data_N(:).group_range];
ground_range = [ray_data_N(:).ground_range];
plot(elevs(idx_goodray), ground_range(idx_goodray), 'b.', 'markersize', 8)
hold on
set(gca, 'fontsize', 14)
grid on


%
% Now do 2D NRT
%
num_range = 600;       
range_inc = 5;
max_range = range_inc * (num_range - 1);
iono_en_grid_2D = zeros(num_ht, num_range);
iono_en_grid_5_2D = zeros(num_ht, num_range);
collision_freq_2D = zeros(num_ht, num_range);
irreg = zeros(4, num_range);
for ii = 1:num_range
  iono_en_grid_2D(:, ii) = en_profile;
  iono_en_grid_5_2D(:, ii) = en_profile_5;
  collision_freq_2D(:, ii) = cf_profile;
end

irregs_flag = 0;
fprintf('2D NRT: generating %d ''no-field'' rays ...', num_elevs);
tic
[ray_data, ray_path_data] = ...
    raytrace_2d(origin_lat, origin_long, elevs, ray_bears(1), freqs, nhops, ...
                tol, irregs_flag, iono_en_grid_2D, iono_en_grid_5_2D, ...
	        collision_freq_2D, ht_start, ht_inc, range_inc, irreg);
toc	    
idx_goodray_2D = find([ray_data(:).ray_label] == 1);
group_range_2D = [ray_data(:).group_range];
ground_range_2D = [ray_data(:).ground_range];

plot(elevs(idx_goodray_2D), ground_range_2D(idx_goodray_2D), '.r', ...
    'markersize', 8)
ylabel('ground range (km)', 'fontsize', 14)
xlabel('elevation (degrees)', 'fontsize', 14)
lh = legend(gca, '3D NRT', '2D NRT');
set(lh, 'fontsize', 14)
hold off

fig1_pos = get(gcf, 'position');
fig2_pos = fig1_pos;
fig1_pos(1) = fig1_pos(1) - 300;
set(gcf, 'position', fig1_pos)

fig2_pos(1) = fig2_pos(1) + 300;


figure(2)
set(gcf, 'position', fig2_pos)
idx = find([ray_data(:).ray_label] == 1 & [ray_data_N(:).ray_label] == 1);  
group_diff = [ray_data_N(:).group_range] - [ray_data(:).group_range];
ground_diff = [ray_data_N(:).ground_range] - [ray_data(:).ground_range];
plot(elevs(idx), ground_diff(idx) * 1000, 'b.')
hold on
plot(elevs(idx), group_diff(idx) * 1000, 'r.')
hold off
set(gca, 'ylim', [-500 500], 'fontsize', 14)
ylabel('3D NRT - 2D NRT range (m)', 'fontsize', 14)
xlabel('elevation (degrees)', 'fontsize', 14)
lh = legend(gca, 'ground range', 'group range');
set(lh, 'fontsize', 14)
grid on
