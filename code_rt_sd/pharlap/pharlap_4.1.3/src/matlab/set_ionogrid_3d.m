%
% Name :
%   set_ionogrid_3d.m adapted from ray_test_3d.m
%
% Purpose :
%   Example of using raytrace_3d for a fan of rays. 
%
% Calling sequence :
%   ray_test_3d
%
% Inputs :
%   None
%
% Outputs :
%   None
%
% Modification History:
%   V1.0  M.A. Cervera  07/12/2009
%     Initial version.
%
%   V1.1  M.A. Cervera  12/05/2009
%     Uses 'parfor' to parallelize the computation if the parallel computing 
%     tool box is available
%
%   V1.3  M.A. Cervera  19/05/2011
%     More efficient handling of ionospheric  and geomagnetic grids grids in
%     call to raytrace_3d  
%
%   V2.0 M.A. Cervera  03/05/2016
%     Modified to make use of multi-threaded raytrace_3d. IRI2016 is now used
%     to generate the ionosphere.
%

%
% setup general stuff
%

% set up ionospheric grid in latitude, longitude and height for uneclipsed
% and eclipsed ionosphere
% 
UT = [2016 8 21 18 30];
UT_array = [2016 8 21 18 30];            % UT - year, month, day, hour, minute
speed_of_light = 2.99792458e8;
R12 = 100;
kp = 3;                         % kp for irregularities model (3 = quiet)
doppler_flag = 1;               % interested in Doppler shift

% range of ray frequencies (MHz)
freq_start= 2;              
freq_inc=2;
freq_end=8;

% input parameters for fixed azimuth and different elevations used by
% plot_elev_3d.m
elevs_a = [45:9:90];               % initial elevation of rays
fixed_azim = 0;                    % azimuth angle
ray_bears_a = fixed_azim * ones(size(elevs_a));
num_elevs = length(elevs_a);

% input parameters for fixed elevation and varyig azimuths used by
% plot_azim_3d.m
ray_bears = [-70:10:90];           % initial bearing of rays
fixed_elev = 60;                   % elevation
elevs = fixed_elev * ones(size(ray_bears));
num_raybears = length(ray_bears);

% Origin point of rays: Fort Hays (Here)
origin_lat = 40.8;              % latitude of the start point of rays
origin_long = -98.336;          % longitude of the start point of rays
origin_ht = 0.0;                % altitude of the start point of rays

%
% generate ionospheric, geomagnetic and irregularity grids
%
ht_start = 100;          % start height for ionospheric grid (km)
ht_inc = 3;              % height increment (km)
num_ht = 200;           
lat_start = 20;          % start latitude at index 1   
lat_inc = 1;             % latitude increment (deg)
num_lat = 100.0;
lon_start= -125;         % start longitude at index 1
lon_inc = 1;             % longitude increment (deg)
num_lon = 100.0;
iono_grid_parms = [lat_start, lat_inc, num_lat, lon_start, lon_inc, num_lon, ...
      ht_start, ht_inc, num_ht, ];

B_ht_start = ht_start;          % start height for geomagnetic grid (km)
B_ht_inc = 10;                  % height increment (km)
B_num_ht = ceil(num_ht .* ht_inc ./ B_ht_inc) + 1;
B_lat_start = lat_start;
B_lat_inc = 1;
B_num_lat = ceil(num_lat .* lat_inc ./ B_lat_inc) + 1;
B_lon_start = lon_start;
B_lon_inc = 1;
B_num_lon = ceil(num_lon .* lon_inc ./ B_lon_inc) + 1; 
geomag_grid_parms = [B_lat_start, B_lat_inc, B_num_lat, B_lon_start, ...
      B_lon_inc, B_num_lon, B_ht_start, B_ht_inc, B_num_ht];
nhops = 1;                  % number of hops
tol = [1e-7 0.01 25];       % ODE solver tolerance and min max stepsizes

cnt =1;
for i=1:size(UT_array,1)
    UT=UT_array(i,:);
tic
fprintf('Generating ionospheric grids... ')
% Grid parameters for uneclipsed, normal ionosphere case
[orig_iono_pf_grid, orig_iono_pf_grid_5, collision_freq, Bx, By, Bz] = ...
    gen_iono_grid_3d_new(UT, R12, iono_grid_parms, ...
                    geomag_grid_parms, doppler_flag, 'iri2016');
toc

% convert plasma frequency grid to  electron density in electrons/cm^3 for
% normal ionosphere
orig_iono_en_grid = orig_iono_pf_grid.^2 / 80.6164e-6;
orig_iono_en_grid_5 = orig_iono_pf_grid_5.^2 / 80.6164e-6;

% Grid parameters for eclipsed ionosphere case, multiply the plasma
% frequency grid parameter with the attenuation factor (lat_lon_new1(:,3))
% at that latitude and longitude
inc=1;
%end_idx = length (lat_lon_new1);
%end_idx = length (lat_lon_new);
end_idx = length (lat_lon_attn);
iono_pf_grid=orig_iono_pf_grid;
iono_pf_grid_5 = orig_iono_pf_grid_5;
% The grid is plama frequency value for a latitude looped over all
% longitudes for all latitudes

% compute the corresponding indexes for iono_pf_grid to apply attenuation
% factor
lat_start_idx = fix((lat_lon_attn(1,1) - lat_start)/lat_inc) + 1; 
lat_end_idx = fix(lat_lon_attn(end_idx,1)-lat_start) + 1;
lon_start_idx = fix((lat_lon_attn(1,2) - lon_start)/lon_inc) + 1; 
lond_end_idx = fix(lat_lon_attn(end_idx,2)-lon_start) + 1;

for lat_i = lat_start_idx:1:lat_end_idx-1 
    for lon_i= lon_start_idx:1:lond_end_idx-1 
        iono_pf_grid(lat_i,lon_i,:) = iono_pf_grid(lat_i,lon_i,:)*(lat_lon_attn(inc,3))^0.5;
        inc=inc+1;
    end
    inc=inc+1;
end

% convert plasma frequency grid to  electron density in electrons/cm^3 for
% eclipsed ionosphere
iono_en_grid = iono_pf_grid.^2 / 80.6164e-6;
iono_en_grid_5 = iono_pf_grid_5.^2 / 80.6164e-6;

% Generate the rays for the case where the magnetic field is ignored 
OX_mode = 0;
% run('NoEclipse_3d_azim_new.m');
% run('Eclipse_3d_azim_new.m');
% script for comparing raytracing normal and eclipsed for fixed azimuth, varing elevations 
% run('plot_elev_3d.m');
run('plot_azim_3d.m');
end


