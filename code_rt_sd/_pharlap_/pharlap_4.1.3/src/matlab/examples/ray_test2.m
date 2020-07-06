%
% Name :
%   ray_test2.m
%
% Purpose :
%   Example of using raytrace_2d for a single ray. Non-default input ray state
%   vector examples are also included. IRI2016 used to generate the ionosphere.
%
% Calling sequence :
%   ray_test2
%
% Inputs :
%   None
%
% Outputs :
%   None
%
% Author:
%   V1.0  M.A. Cervera  16/06/2006
%
%   V1.1  M.A. Cervera  28/09/2006  
%     Minor update against raytrace_2d
%
%   V2.0  M.A. Cervera  15/06/2007  
%     Updated for Pharlap V2. Extra raytracing calls to exemplify the new ray 
%     state vector features of raytrace_2d.
%
%   V2.1 M.A. Cervera  12/03/2008
%      Minor update against raytrace_2d. 
%
%   V2.2 M.A. Cervera  01/05/2008
%      Renamed to ray_test2 and modified to work with updated raytrace_2d
%      (for pre-release version of PHaRLAP 3.0)
%
%   V2.3 M.A. Cervera  15/12/2008
%      Now uses IRI2007 to generate the ionosphere
%
%   V2.4 M.A. Cervera  19/05/2011
%      More efficient handling of ionospheric grids in call to raytrace_2d 
%
%   V2.5  M.A. Cervera  02/05/2016
%      Updated to use IRI2016
%
%   V2.6  M.A. Cervera  20/05/2016
%      Updated to use multi-threaded raytrace_2d
%

%
% setup general stuff
%
UT = [2000 9 21 0 0];        % UT - year, month, day, hour, minute
speed_of_light = 2.99792458e8;
R12 = 100;
elev = 10;                   % ray elevation
ray_bear = 329;              % bearing of ray 
origin_lat = -20.0;          % latitude of the start point of ray
origin_long = 130.0;         % longitude of the start point of ray
doppler_flag = 0;            % not interested in Doppler shift
irregs_flag = 0;             % no irregularities
kp = 0;                      % kp not used as irregs_flag = 0. Set it to a 
                             % dummy value 

%
% generate ionospheric, geomagnetic and irregularity grids
%
max_range = 10000;      % maximum range for sampling the ionosphere (km)
num_range = 201;        % number of ranges (must be < 2001)
range_inc = max_range ./ (num_range - 1);  % range cell size (km)

start_height = 60;      % start height for ionospheric grid (km)
height_inc = 2;         % height increment (km)
num_heights = 100;      % number of  heights (must be < 2001)

clear iri_options
iri_options.Ne_B0B1_model = 'Bil-2000'; % this is a non-standard setting for 
                                        % IRI but is used as an example

tic
[iono_pf_grid, iono_pf_grid_5, collision_freq, irreg] = ...
   gen_iono_grid_2d(origin_lat, origin_long, R12, UT, ray_bear, ...
                    max_range, num_range, range_inc, start_height, ...
		    height_inc, num_heights, kp, doppler_flag, 'iri2016', ...
		    iri_options);
toc

% convert plasma frequency grid to  electron density in electrons/cm^3
iono_en_grid = iono_pf_grid.^2 / 80.6164e-6;
iono_en_grid_5 = iono_pf_grid_5.^2 / 80.6164e-6;

%
% call raytrace
%
freq = 20.0;         % frequency (MHz)
nhops = 1;           % number of hops
elev = 10;
tol = [1e-7 .001 100]; % ODE tolerance and min/max step sizes
tic
tol = 1e-7;

[ray_data, ray_path_data, ray_state_vec] = ...
      raytrace_2d(origin_lat, origin_long, elev, ray_bear, freq, nhops, tol, ...
             irregs_flag, iono_en_grid, iono_en_grid_5, collision_freq, ...
	     start_height, height_inc, range_inc, irreg);
	 	 
toc    
% plot the ray	 
plot(ray_path_data.ground_range, ray_path_data.height, 'b*')
set(gca,'xlim', [0 2000])
xlabel('Ground Range (km)')
ylabel('Height (km)')

%
% new raytrace using a non-default input ray state vector - use that from a
% point along the previous raytrace (90 km in altittude) and see if we can
% replicate the previous raytrace 
%
idx = find(ray_path_data.height > 90); idx = idx(1);
ray_state_vec_in2.r = ray_state_vec.r(idx);
ray_state_vec_in2.Q = ray_state_vec.Q(idx);
ray_state_vec_in2.theta = ray_state_vec.theta(idx);
ray_state_vec_in2.delta_r = ray_state_vec.delta_r(idx);
ray_state_vec_in2.delta_Q = ray_state_vec.delta_Q(idx);
ray_state_vec_in2.deviative_absorption = ...
                ray_state_vec.deviative_absorption(idx);
ray_state_vec_in2.phase_path = ray_state_vec.phase_path(idx);
ray_state_vec_in2.group_path = ray_state_vec.group_path(idx);
ray_state_vec_in2.group_step_size = ray_state_vec.group_step_size(idx);

[ray_data2, ray_path_data2] = ...
    raytrace_2d(origin_lat, origin_long, elev, ray_bear, freq, nhops, ...
             tol, irregs_flag, ray_state_vec_in2);

% overplot the ray	 
hold on
plot(ray_path_data2.ground_range, ray_path_data2.height,'*r')
hold off



%
% 3rd raytrace using a non-default input ray state vector - use that from a
% point along the previous raytrace (90 km in altittude) and see if we can
% reflect the previous raytrace 
%
ray_state_vec_in3 = ray_state_vec_in2;
ray_state_vec_in3.Q = -ray_state_vec_in3.Q;
[ray_data3, ray_path_data3] = ...
    raytrace_2d(origin_lat, origin_long, elev, ray_bear, freq, nhops, ...
             tol, irregs_flag, ray_state_vec_in3);

% overplot the ray	 
hold on
plot(ray_path_data3.ground_range, ray_path_data3.height,'*g')
plot([300 600], [90 90], 'k')
text(100, 90, 'Es Layer')
hold off


%
% display the legend
%
legend('Original ray', 'Ray unaffected by Es Layer', ...
       'Ray reflected by Es Layer')