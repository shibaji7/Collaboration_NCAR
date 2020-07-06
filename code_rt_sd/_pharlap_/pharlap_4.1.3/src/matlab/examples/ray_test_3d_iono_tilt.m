%
% Name :
%   ray_test_3d_iono_tilt.m
%
% Purpose :
%   Example of 3D raytracing through a tilted ionosphere.
%
% Calling sequence :
%   ray_test_3d_iono_tilt
%
% Inputs :
%   None.
%
% Outputs :
%   None.
%
% Author:
%   V1.0  M.A. Cervera  30/09/2010
%
%   V1.2  M.A. Cervera  19/05/2011
%      More efficient handling of ionospheric and geomagnetic grids in call
%      to raytrace_3d 
%
%   V2.0  M.A. Cervera  03/05/2016
%      Updated to use the multi threaded version of raytrace_3d. IRI2016 is
%      now used to generate the ionosphere.
%


%
% setup general stuff
%
UT = [2008 2 26 6 0];              % UT - year, month, day, hour, minute
R12 = 70;                          % yearly smoothed sunspot number
pfsq_conv = 80.6163849431291e-12;  % multiplicative factor to convert electron
                                   % dens (m^3) to plasma freq. squared (MHz^2)
tx_lat = -35;                      % tx. location
tx_lon = 138;
tx_ht = 0.0;
ray_bearing = 0.0;                 % ray bearing


%
% Generate the 3D ionospheric, geomagnetic and irregularity grids
%
ht_start = 60;          % start height for ionospheric grid (km)
ht_inc = 2.5;           % height increment (km)
num_ht = 200;           % number of  heights (must be <= 401)
lat_start = -36.0;
lat_inc = 0.3;
num_lat = 20;           % number of latitudes (must be <= 701)
lon_start= 135;
lon_inc = 0.3;
num_lon = 20;           % number of longitudes (must be <= 701)
iono_grid_parms = [lat_start, lat_inc, num_lat, lon_start, lon_inc, num_lon, ...
      ht_start, ht_inc, num_ht, ];

B_ht_start = ht_start;          % start height for geomagnetic grid (km)
B_ht_inc = 10;                  % height increment (km)
B_num_ht = ceil(num_ht .* ht_inc ./ B_ht_inc) + 1;
B_lat_start = lat_start;
B_lat_inc = 0.5;
B_num_lat = ceil(num_lat .* lat_inc ./ B_lat_inc) + 1;
B_lon_start = lon_start;
B_lon_inc = 0.5;
B_num_lon = ceil(num_lon .* lon_inc ./ B_lon_inc) + 1; 
geomag_grid_parms = [B_lat_start, B_lat_inc, B_num_lat, B_lon_start, ...
      B_lon_inc, B_num_lon, B_ht_start, B_ht_inc, B_num_ht];

doppler_flag = 0;

fprintf('Generating 3D ionospheric and geomag grids... ')
tic
[iono_pf_grid, iono_pf_grid_5, collision_freq, Bx, By, Bz] = ...
    gen_iono_grid_3d(UT, R12, iono_grid_parms, geomag_grid_parms, ...
                     doppler_flag, 'iri2016');

% make the ionosphere "spherically symmetric" - use same ionospheric profile
% at all lat, lons
[iono, iono_extra] = iri2012(tx_lat, tx_lon, R12, UT);	
foE  = sqrt(iono_extra(5) .* pfsq_conv);    
hmE  = iono_extra(6); 
if (iono_extra(3) == -1)
  foF1 = 0;
  hmF1 = 0;
else
  foF1 = sqrt(iono_extra(3) .* pfsq_conv); 
  hmF1 = iono_extra(4); 
end
foF2 = sqrt(iono_extra(1) .* pfsq_conv);    
hmF2 = iono_extra(2);
% layer thicknesses - IRI does not return these values - use FIRIC model
ymE  = 20.0; 
ymF1 = hmF1 ./ 4;
sol_zen_ang = solar_za(tx_lat, tx_lon, UT);
SZA_term = 0.15 .* exp(-(sol_zen_ang./105).^4) + 0.25;	
ymF2 = SZA_term .* hmF2;

ht_max = ht_start + (num_ht - 1) * ht_inc;
height_arr = [ht_start : ht_inc : ht_max];
for lonidx = 1:num_lon
  foF2_mod = foF2;
  iono_profile = chapman(foE, hmE, ymE, foF1, hmF1, ymF1, foF2_mod, hmF2, ...
                         ymF2, height_arr);
  for latidx = 1:num_lat
    iono_pf_grid(latidx, lonidx, :) = iono_profile;
    iono_pf_grid_5(latidx, lonidx, :) = iono_profile;
  end
end

toc

% convert plasma frequency grid to  electron density in electrons/cm^3
iono_en_grid = iono_pf_grid.^2 / 80.6164e-6;
iono_en_grid_5 = iono_pf_grid_5.^2 / 80.6164e-6;


%
% Now for the 3D raytracing through the "spherically symmetric" ionosphere
%
fprintf('Do the 3D raytracing... \n')
group_range_O = [];
freqs_O = [];
OX_mode = 1;
tol = 1e-7;
nhops = 1;
els = [50 : 2 : 80 ]; 
azs = 30*[-1.0 : 0.5 : 1.0];
num_elevs = length(els);
num_azims = length(azs);
num_rays = num_elevs.*num_azims;
elevs = repmat(els, 1, num_azims);
azims = [];
for ii=1:num_azims
  azims = [azims ones(1, num_elevs).*azs(ii)];
end
freqs = ones(1, num_rays).*5;
  
lat = [];
lon = [];
grp_rng = [];
el = [];
az = [];
ray_labels = [];

[ray_data, ray_path_data, ray_state_vec] = ...
	 raytrace_3d(tx_lat, tx_lon, tx_ht, elevs, azims, freqs, ...
	     OX_mode, nhops, tol, iono_en_grid, iono_en_grid_5, ...
	     collision_freq, iono_grid_parms, Bx, By, Bz, ...
	     geomag_grid_parms);

clear raytrace_3d    % don't need this ionosphere again so clear it

figure(2)	 
az_0_idx = find(azims == 0);
for ii = 1:length(az_0_idx)
  idx = az_0_idx(ii);
  plot3(ray_path_data(idx).lon, ray_path_data(idx).lat, ...
        ray_path_data(idx).height, 'b.')
  hold on
end
grid on
set(gca, 'linewidth', 2, 'fontsize', 14)
xlabel('Longitude  ')
ylabel('Latitude  ')
zlabel('Height  ')

figure(1)
lats = zeros(1, num_rays) .* NaN;
lons = zeros(1, num_rays) .* NaN;
for ii = 1:num_rays
  if (ray_data(ii).ray_label == 1)
    lats(ii)= ray_data(ii).lat;
    lons(ii) = ray_data(ii).lon;
  end
end
plot(lons, lats, 'o', 'markerfacecolor', 'b', 'markeredgecolor', 'b');
set(gca, 'fontsize', 14, 'linewidth', 2, 'ylim', [-35.2 -31], ...
    'xlim', [135.5 140.5])
xlabel('Longitude', 'fontsize', 14)
ylabel('Latitude', 'fontsize', 14)
hold on
plot(tx_lon, tx_lat, 'ko', 'markerfacecolor', 'k');



%
% Now introduce tilt in the ionosphere - increase foF2 as function of longitude
% 0.03MHz increase per degree of long.
for lonidx = 1:num_lon
  foF2_mod = foF2 + 2*lonidx/num_lon;
  iono_profile = chapman(foE, hmE, ymE, foF1, hmF1, ymF1, foF2_mod, hmF2, ...
                         ymF2, height_arr);
  for latidx = 1:num_lat
    iono_pf_grid(latidx, lonidx, :) = iono_profile;
    iono_pf_grid_5(latidx, lonidx, :) = iono_profile;
  end
end

% convert plasma frequency grid to  electron density in electrons/cm^3
iono_en_grid = iono_pf_grid.^2 / 80.6164e-6;
iono_en_grid_5 = iono_pf_grid_5.^2 / 80.6164e-6;



%
% Now for the 3D raytracing with an ionospheric tilt
%
[ray_data, ray_path_data, ray_state_vec] = ...
	 raytrace_3d(tx_lat, tx_lon, tx_ht, elevs, azims, freqs, ...
	     OX_mode, nhops, tol, iono_en_grid, iono_en_grid_5, ...
	     collision_freq, iono_grid_parms, Bx, By, Bz, ...
	     geomag_grid_parms);

clear raytrace_3d    % don't need this ionosphere again so clear it
	 
figure(2)	 
az_0_idx = find(azims == 0);
for ii = 1:length(az_0_idx)
  idx = az_0_idx(ii);
  plot3(ray_path_data(idx).lon, ray_path_data(idx).lat, ...
        ray_path_data(idx).height, 'g.')
  hold on
end
grid on
set(gca, 'linewidth', 2, 'fontsize', 14)
xlabel('Longitude  ')
ylabel('Latitude  ')
zlabel('Height  ')

figure(1)
lats = zeros(1, num_rays) .* NaN;
lons = zeros(1, num_rays) .* NaN;
for ii = 1:num_rays
  if (ray_data(ii).ray_label == 1)
    lats(ii)= ray_data(ii).lat;
    lons(ii) = ray_data(ii).lon;
  end
end
plot(lons, lats, 'o', 'markerfacecolor', 'g', 'markeredgecolor', 'g');
set(gca, 'fontsize', 14, 'linewidth', 2, 'ylim', [-35.2 -31], ...
    'xlim', [135.5 140.5])
xlabel('Longitude', 'fontsize', 14)
ylabel('Latitude', 'fontsize', 14)
hold on
plot(tx_lon, tx_lat, 'ko', 'markerfacecolor', 'k');



%
% Now increase the tilt in ionosphere to 0.075MHz per degree in longitude
%
for lonidx = 1:num_lon
  foF2_mod = foF2 + 5*lonidx/num_lon;
  iono_profile = chapman(foE, hmE, ymE, foF1, hmF1, ymF1, foF2_mod, hmF2, ...
                         ymF2, height_arr);
  for latidx = 1:num_lat
    iono_pf_grid(latidx, lonidx, :) = iono_profile;
    iono_pf_grid_5(latidx, lonidx, :) = iono_profile;
  end
end

% convert plasma frequency grid to  electron density in electrons/cm^3
iono_en_grid = iono_pf_grid.^2 / 80.6164e-6;
iono_en_grid_5 = iono_pf_grid_5.^2 / 80.6164e-6;



%
% Perfrom the 3D raytracing with the larger ionospheric tilt
%
[ray_data, ray_path_data, ray_state_vec] = ...
	 raytrace_3d(tx_lat, tx_lon, tx_ht, elevs, azims, freqs, ...
	     OX_mode, nhops, tol, iono_en_grid, iono_en_grid_5, ...
	     collision_freq, iono_grid_parms, Bx, By, Bz, ...
	     geomag_grid_parms);

clear raytrace_3d    % don't need this ionosphere again so clear it
	 
figure(2)	 
az_0_idx = find(azims == 0);
for ii = 1:length(az_0_idx)
  idx = az_0_idx(ii);
  plot3(ray_path_data(idx).lon, ray_path_data(idx).lat, ...
        ray_path_data(idx).height, 'r.')
  hold on
end
hold off
grid on
set(gca, 'linewidth', 2, 'fontsize', 14)
xlabel('Longitude  ')
ylabel('Latitude  ')
zlabel('Height  ')

figure(1)
lats = zeros(1, num_rays) .* NaN;
lons = zeros(1, num_rays) .* NaN;
for ii = 1:num_rays
  if (ray_data(ii).ray_label == 1)
    lats(ii)= ray_data(ii).lat;
    lons(ii) = ray_data(ii).lon;
  end
end
plot(lons, lats, 'o', 'markerfacecolor', 'r', 'markeredgecolor', 'r');
set(gca, 'fontsize', 14, 'linewidth', 2, 'ylim', [-35.2 -31], ...
    'xlim', [135.5 140.5])
xlabel('Longitude', 'fontsize', 14)
ylabel('Latitude', 'fontsize', 14)
hold on
plot(tx_lon, tx_lat, 'ko', 'markerfacecolor', 'k');



figure(1)
hold off

figure(2)
hold off
