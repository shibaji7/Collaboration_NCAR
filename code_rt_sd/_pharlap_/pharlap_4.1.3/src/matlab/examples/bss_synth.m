%
% Name :
%   bss_synth.m
%
% Purpose :
%   Synthetic back-scatter ionogram generation using 2d raytracing. Antenna
%   gain patterns and receiver miss-match loss are not included (search for
%   gain_tx, gain_rx and rx_miss_match_loss to see where they should be
%   included). 
%
% Calling sequence :
%   bss_synth
%
% Inputs :
%   
%
% Outputs :
%   
%
% Author:
%   V1.0  M.A. Cervera  16/04/2008
%
%   V1.1  M.A. Cervera  01/05/2008
%      Modified to work with updated raytrace_2d (for pre-release version of
%      PHaRLAP 3.0)
%
%   V1.2 M.A. Cervera  11/09/2009
%      Corrected error with calculation of D region absorption
%
%   V1.3 M.A. Cervera  20/05/2014
%      Improved handling of flux tubes.
%
%   V1.4  M.A. Cervera  02/05/2016
%      Updated to use IRI2016
%
%   V1.5  M.A. Cervera  20/05/2016
%      Updated to use multi-threaded raytrace_2d
%


%
% setup physical constants and parameters
%
radius_earth = 6371;          % Mean radius of the Earth (km)
light_speed = 299792459.0;    % speed of light 

%
% setup general stuff
%
UT = [2001 6 15 9 30];        % UT - year, month, day, hour, minute
R12 = 130;                    % yearly smoothed sunspot number

tx_lat = -25;                 % latitude of the start point of ray 
tx_long = 135;                % longitude of the start point of ray 

tx_pow = 10;                  % transmitter power (Watts)
tx_pow = 10 .* log10(tx_pow); % transmitter power (dBW)

ray_bearing = 310;            % direction (degrees) of the beam

num_elem_rx_array = 15;       % no. of antenna elements in rx. array
rx_arr_space = 6;             % spacing between antenna elements
arr_length = (num_elem_rx_array - 1) .* rx_arr_space;

%
% generate ionospheric, geomagnetic and irregularity grids
%
max_range = 8000;       % maximum range for sampling the ionosphere (km)
num_range = 201;        % number of ranges (must be < 2001)
range_inc = max_range ./ (num_range - 1);  % range cell size (km)
range_depth = range_inc;

start_height = 60;      % start height for ionospheric grid (km)
height_inc = 5;         % height increment (km)
num_heights = 201;      % number of  heights (must be < 2001)

doppler_flag = 0;       % not interested in Doppler shift and spread
kp = 0;                 % kp not used as doppler_flag = 0. Set it to a 
                        % dummy value 

clear iri_options			
iri_options.Ne_B0B1_model = 'Bil-2000';
[iono_pf_grid, iono_pf_grid_5, collision_freq, irreg] = ...
    gen_iono_grid_2d(tx_lat, tx_long, R12, UT, ray_bearing, ...
                     max_range, num_range, range_inc, start_height, ...
		     height_inc, num_heights, kp, doppler_flag, 'iri2016', ...
		     iri_options);

% convert plasma frequency grid to electron density in electrons/cm^3
iono_en_grid = iono_pf_grid.^2 / 80.6164e-6;
iono_en_grid_5 = iono_pf_grid_5.^2 / 80.6164e-6;


%
% set up necessary parameters for the raytacing
%
irregs_flag = 0;     % no irregularities
num_hops = 4;        % number of hops
tol = 1e-7;          % rkf tolerance

freq_step = 0.5;
freq_min = 2;
freq_max = 45;
elev_step = 0.5;
elev_min = 2;
elev_max = 50;
freqs = (freq_min: freq_step: freq_max);
num_freqs = length(freqs);
elevs = (elev_min : elev_step : elev_max);
num_elevs = length(elevs);

group_step = 50;
group_max = 10000;
groups = (0: group_step : group_max);
num_groups = length(groups);

range_res = group_step;

power_total = zeros(num_groups, num_freqs) + 1e-200;
power_dominant = zeros(num_groups, num_freqs);
phase_path_dom = zeros(num_groups, num_freqs);
elevation_dom = zeros(num_groups, num_freqs);
num_hops_dom = zeros(num_groups, num_freqs);

% If it doesn't already exist, set up the figure for plotting the ionogram 
% as it is being synthesised
if exist('bss_fig_handle', 'var')
  if ~ishandle(bss_fig_handle)
    initialise_figure = true; 
  else
    initialise_figure = false;
  end
else
  initialise_figure = true;
end

if initialise_figure
  bss_fig_handle = figure('Name', 'Synthetic BSS Ionogram', ...
    'NumberTitle', 'off', 'Position', [150 500 750 550]);
  xrange = [freq_min freq_max];
  yrange = [0 group_max];
  bss_image_handle = imagesc(xrange, yrange, power_total, [-300 -100]);
  set(gca, 'ydir', 'normal', 'Fontsize', 14)
  ylabel('Group Range (km)', 'Fontsize', 14)
  xlabel('Frequency (MHz)', 'Fontsize',  14)
  cb_handle = colorbar('Fontsize', 14);
  set(get(cb_handle,'ylabel'), 'string', 'Power (dB)', 'Fontsize', 14, ...
    'verticalalignment','top')
  
end



% first call to raytrace so pass in the ionospheric and geomagnetic grids 
[ray_data, ray_path_data] = ...
    raytrace_2d(tx_lat, tx_long, elevs(1), ray_bearing, freq_min, ...
	       num_hops, tol, irregs_flag, iono_en_grid, iono_en_grid_5, ...
	       collision_freq, start_height, height_inc, range_inc, irreg);


%
% loop over ray frequency
%
for freq = freq_min:freq_step:freq_max

  freq_bin = (freq - freq_min) ./ freq_step + 1; 

  latitude_hop = zeros(num_elevs, num_hops) * nan;
  longitude_hop = zeros(num_elevs, num_hops) * nan;
  gnd_range = zeros(num_elevs, num_hops) * nan;
  grp_range = zeros(num_elevs, num_hops) * nan;
  elev_init = zeros(num_elevs, num_hops) * nan;
  elev_fin = zeros(num_elevs, num_hops) * nan;
  elevation_tx = zeros(num_elevs, num_hops) * nan;
  effective_range = zeros(num_elevs, num_hops) * nan;
  dev_absorption = zeros(num_elevs, num_hops) * nan;
  nondev_absorption = zeros(num_elevs, num_hops) * nan;
  forward_scatt_loss = zeros(num_elevs, num_hops) * nan;
  phase_path = zeros(num_elevs, num_hops) * nan;
  labels = zeros(num_elevs, num_hops) * nan;
  num_hops_done = zeros(num_elevs, 1);
  
  % calculate the azimuthal beam width
  wave_length = light_speed ./ (freq *1e6);             % wave length in m
  beam_az_width = wave_length ./ arr_length;            % azimuthal beam width

  % call 2D raytrace - multi-thread over elevation
  freqs = freq .* ones(size(elevs));
  [ray_data, ray_path_data] = ...
	raytrace_2d(tx_lat, tx_long, elevs, ray_bearing, freqs, ...
	            num_hops, tol, irregs_flag);
 
  % loop over ray elevation and calculate various quantites
  for elev_idx = 1:num_elevs 
    nhops_done = ray_data(elev_idx).nhops_attempted;
          
    % obtain the ground and group range, initial and final elevation,
    % effective range, deviative absorption, phase path, O-X mode correction,
    % number of hops done and hop label.
    latitude_hop(elev_idx, 1:nhops_done) = ray_data(elev_idx).lat;
    longitude_hop(elev_idx, 1:nhops_done) = ray_data(elev_idx).lon;
    gnd_range(elev_idx, 1:nhops_done) = ray_data(elev_idx).ground_range;
    grp_range(elev_idx, 1:nhops_done) = ray_data(elev_idx).group_range;
    elev_init(elev_idx, 1:nhops_done) = ray_data(elev_idx).initial_elev;
    elev_fin(elev_idx, 1:nhops_done) = ray_data(elev_idx).final_elev;
    elevation_tx(elev_idx, 1:nhops_done) = elevs(elev_idx);
    effective_range(elev_idx, 1:nhops_done) = ...
	                       ray_data(elev_idx).effective_range; 
    dev_absorption(elev_idx, 1:nhops_done) = ...
	                       ray_data(elev_idx).deviative_absorption;
    phase_path(elev_idx, 1:nhops_done) = ray_data(elev_idx).phase_path;
    num_hops_done(elev_idx) = ray_data(elev_idx).nhops_attempted;
    labels(elev_idx, 1:nhops_done) = ray_data(elev_idx).ray_label;
 
    % calculate the cumulative-hop forward scattering loss and ionospheric
    % non-deviative absorption loss (O mode only)
    O_mode = 1;
    fsloss_hop_cumu = 0;
    absorp_cumu = 0;
    lat_hop_start = tx_lat;
    lon_hop_start = tx_long;
    elev_hop_start = elevs(elev_idx);
    forward_scatt_loss(elev_idx, 1) = 0;
    for hop_idx = 1:nhops_done-1
      lat_hop_end = ray_data(elev_idx).lat(hop_idx);
      lon_hop_end = ray_data(elev_idx).lon(hop_idx);   
      elev_hop_end = elev_fin(elev_idx, hop_idx);
      
      % forward scattering loss
      fsloss_hop_cumu = fsloss_hop_cumu + ...
	  ground_fs_loss(lat_hop_end, lon_hop_end, elev_hop_end, freq);
      forward_scatt_loss(elev_idx, hop_idx+1) = fsloss_hop_cumu;
      
      % absorption loss 
      apogee_gndr = ray_data(elev_idx).gnd_rng_to_apogee(hop_idx);
      [lat_hop_mp, lon_hop_mp] = raz2latlon(apogee_gndr*1000, ray_bearing, ...
	    lat_hop_start, lon_hop_start);
      absorp_thishop = abso_bg(lat_hop_mp, lon_hop_mp, elev_hop_start, ...
	    freq, UT, R12, O_mode);
      absorp_cumu = absorp_cumu + absorp_thishop;
      nondev_absorption(elev_idx, hop_idx) = absorp_cumu;
      
      lat_hop_start = lat_hop_end;
      lon_hop_start = lon_hop_end; 
      elev_hop_start = elev_hop_end;
    end
    
    % add on the absorption loss for the last hop
    hop_idx = hop_idx + 1;
    if labels(elev_idx, hop_idx) == 1
      apogee_gndr = ray_data(elev_idx).gnd_rng_to_apogee(hop_idx);
      [lat_hop_mp, lon_hop_mp] = raz2latlon(apogee_gndr*1000, ray_bearing, ...
  	  lat_hop_start, lon_hop_start);
      absorp_lasthop = abso_bg(lat_hop_mp, lon_hop_mp, elev_hop_start, ...
  	                     freq, UT, R12, O_mode);
      absorp_cumu = absorp_cumu + absorp_lasthop;
      nondev_absorption(elev_idx, hop_idx) = absorp_cumu;
    end
    
  end   % for elev_idx = 1:num_elevs
 
  % add the non-deviative and dviative absorptions to give the total
  % ionospheric absorption
  iono_absorption = nondev_absorption + dev_absorption;
  
  
  %
  % calculate rate of change of group-range wrt elevation and rate of change
  % of ground-range wrt group-range
  %
  dgrp_dels = zeros(size(grp_range)) .* NaN;
  dgnd_dgrp = zeros(size(grp_range)) .* NaN;
  for ii = 1:num_hops
    good_ray = find(labels(:,ii) > -2);  % these rays have not penetrated
    if length(good_ray) > 2 
      dgrp_dels(good_ray, ii) = deriv(grp_range(good_ray, ii)', ...
        elev_init(good_ray, ii)'); 
      dgnd_dgrp(good_ray, ii) = deriv(gnd_range(good_ray, ii)', ...
	  grp_range(good_ray, ii)');      
    end
  end  
  
  
  %
  % Find the the rays which return to ground. Ignore the other rays as they 
  % don't propogate back to the receiver. Remove rays whose rate of change of
  % range with elevation is too large as this indicates we are too far into a 
  % cusp region to be reliable.
  %
  idx_raygnd = find(labels == 1 & grp_range < group_max & abs(dgrp_dels) < 400);
  idx_reflect = find(labels == 0 & grp_range < group_max);
  idx_fai = find(labels == -1 & grp_range < group_max);

  lat_hop = latitude_hop(idx_raygnd);
  lon_hop = longitude_hop(idx_raygnd);
  gnd_rng = gnd_range(idx_raygnd);
  grp_rng = grp_range(idx_raygnd);
  dgrp_dels = dgrp_dels(idx_raygnd);
  dgnd_dgrp = dgnd_dgrp(idx_raygnd);
  elev_tx = elevation_tx(idx_raygnd);
  phase = phase_path(idx_raygnd);
  fs_loss = forward_scatt_loss(idx_raygnd);
  eff_range = effective_range(idx_raygnd);
  iono_absorp = iono_absorption(idx_raygnd);
  hop_num_array = repmat((1:num_hops),num_elevs,1);
  hop_num = hop_num_array(idx_raygnd);
  
  
  %
  % Loop over all the out-bound rays - find matching in-bound rays (same
  % ground range) and use the radar equation to calculate backscattered
  % power. Bin according to group-range and frequency.
  %
  for ray_out = 1:length(gnd_rng)
    this_gnd_rng = gnd_rng(ray_out);
     
    % various out-bound quantities
    lat_hop_out = lat_hop(ray_out);
    lon_hop_out = (ray_out);
    grp_out = grp_rng(ray_out);
    elev_out = elev_tx(ray_out);
    fs_loss_out = fs_loss(ray_out);
    eff_range_out = eff_range(ray_out);
    iono_absorp_out = iono_absorp(ray_out);
    dgnd_dgrp_out = dgnd_dgrp(ray_out);
    dgrp_dels_out = dgrp_dels(ray_out);
    
    % for this out-bound ray find the matching in-bound rays ie. those which 
    % have the same ground range - we need to find the bracketing rays but 
    % ignore those whose rate of change of range with elevation is too large 
    % as this indicates we are too far into a cusp region to be reliable
    for ray_in=1:length(gnd_rng)-1  

      if ( ( (gnd_rng(ray_in) >= this_gnd_rng && ...
              gnd_rng(ray_in+1) <= this_gnd_rng) || ...
             (gnd_rng(ray_in) <= this_gnd_rng && ...   
              gnd_rng(ray_in+1) >= this_gnd_rng)) && ...
	   (hop_num(ray_in) == hop_num(ray_in+1)) )  
      
        idx = [ray_in ray_in+1];
        grp_in = ...
          interp1(gnd_rng(idx), grp_rng(idx), this_gnd_rng);
        elev_in = ...
          interp1(gnd_rng(idx), elev_tx(idx), this_gnd_rng);
        phase_in = ...
          interp1(gnd_rng(idx), phase(idx), this_gnd_rng);
        iono_absorp_in = ...
          interp1(gnd_rng(idx), iono_absorp(idx), this_gnd_rng);
        eff_range_in = ...
          interp1(gnd_rng(idx), eff_range(idx), this_gnd_rng);
        fs_loss_in = ...
          interp1(gnd_rng(idx), fs_loss(idx), this_gnd_rng);
        dgrp_dels_in = ...
	  interp1(gnd_rng(idx), dgrp_dels(idx), this_gnd_rng);
        num_hops_in = hop_num(ray_in);
     
	group_bin = round((grp_out + grp_in) ./ 2 ./ range_res) + 1;

	% calculate the ground backscatter loss
	bs_loss = ground_bs_loss(lat_hop_out, lon_hop_out);

	% Calculate the area of ground illuminated by the ray - note that
	% dgnd_dgrp -> inf at the "time-caustic" which means that geometrical
	% optics breaks down and the equation used for area is no longer
	% valid. Apply a limit to dgnd_dgrp to account for this. See Ong,
	% Dyson and Bennet, RS, 1173-1186, 1998 for a better solution at the
	% time-caustic. 
	area = radius_earth .* sin(this_gnd_rng ./ radius_earth) .* ...
	       range_depth .* min(abs(dgnd_dgrp_out), 2) .* beam_az_width;
	area = area * 1e6;     % convert units from km^2 to m^2

	% need antenna gains and receiver miss-match loss
	gain_tx = 1;
	gain_rx = 1;
	rx_miss_match_loss = 0;   % dB

	% basic two-way radar equation for ground backscatter of flux-tube
	pow = tx_pow * gain_tx .* gain_rx .* (wave_length.^2 ./ (4.*pi)) .* ...
	      area ./ (16.*pi.^2 .* eff_range_in.^2  * eff_range_out.^2);
	pow_dB = 10 .* log10(pow) - bs_loss;

	% include forward scattering, ionospheric absorption, receiver
	% miss-match, and cone-effect spreading losses
	ray_power_dB = pow_dB - fs_loss_out - fs_loss_in - ...
	    rx_miss_match_loss - iono_absorp_out - iono_absorp_in;

	% calculate group range of mode and the group-range extent of flux-tube
	% represented by the ray
	group = (grp_out + grp_in) ./ 2;
	del_elev_out = elev_step;
	delta_group = del_elev_out .* abs(dgrp_dels_in + dgrp_dels_out) ./ 2;
	group_start = group - delta_group ./ 2;
	group_start = max([1, group_start]);
	group_end = group + delta_group ./ 2;
	group_end = min([group_end group_max]);

	% determine the group range bins to which the mode contributes energy
	group_start_bin = fix(group_start ./ group_step) + 1;
	group_end_bin = fix(group_end ./ group_step) + 1;
	grp_bin_idx = (group_start_bin : 1 : group_end_bin);
	num_bins = length(grp_bin_idx);

	% determine scaling factor for power in each group bin - will be
	% 1/numbins for each bin except at the ends where it is give by the
	% fractional conrtibution and is < 1/numbins
	bin_scale = ones(size(grp_bin_idx));
	if length(bin_scale) > 1
	  bin_scale(1) = group_start_bin - group_start ./ group_step;
	  bin_scale(num_bins) = group_end  ./ group_step - (group_end_bin - 1);
	  bin_scale = bin_scale ./ sum(bin_scale);
	else
	  bin_scale = 1;
	end
        bin_scale = bin_scale';
	
	% populate the ionogram range-frequency grid
	ray_power = 10.^(ray_power_dB ./ 10);
	power_total(grp_bin_idx, freq_bin) = ...
		  power_total(grp_bin_idx, freq_bin) + ray_power.* bin_scale;
	      
	for ii = 1:length(bin_scale)
          if ray_power.*bin_scale(ii) > power_dominant(grp_bin_idx(ii),freq_bin)
  	    power_dominant(grp_bin_idx(ii), freq_bin) =ray_power.*bin_scale(ii);
  	    phase_path_dom(grp_bin_idx(ii), freq_bin) = phase_in;
  	    elevation_dom(grp_bin_idx(ii), freq_bin) = elev_in;
  	    num_hops_dom(grp_bin_idx(ii), freq_bin) = num_hops_in;
          end
        end
 
      end

    end

  end
  
  % plot the ionogram as it is being synthesised
  pow_tot_dB = 10.*log10(power_total);
  pow_dom_dB = 10.*log10(power_dominant);
  set(bss_image_handle, 'CData', pow_tot_dB)
  figure(bss_fig_handle)
  drawnow
  
end   % of frequency loop


%
% plot the dominant mode elevation
%
if exist('dompow_fig_handle', 'var')
  if ~ishandle(dompow_fig_handle)
    initialise_figure = true; 
  else
    initialise_figure = false;
  end
else
  initialise_figure = true;
end

dompow_fig_handle = figure('Name', ...
  'Synthetic BSS Ionogram: Dominant mode power', ...
  'NumberTitle', 'off', 'Position', [150 200 750 550]);
xrange = [freq_min freq_max];
yrange = [0 group_max];
dom_image_handle = imagesc(xrange, yrange, pow_dom_dB, [-300 -100]);
set(gca, 'ydir', 'normal', 'Fontsize', 14)
ylabel('Group Range (km)', 'Fontsize', 14)
xlabel('Frequency (MHz)', 'Fontsize',  14)
cb_handle = colorbar('Fontsize', 14);
set(get(cb_handle,'ylabel'), 'string', 'Power (dB)', 'Fontsize', 14, ...
  'verticalalignment','top')


%
% plot the dominant mode elevation
%
domelev_fig_handle = figure('Name', ...
  'Synthetic BSS Ionogram: Dominant mode elevation', ...
  'NumberTitle', 'off', 'Position', [950 500 750 550]);
xrange = [freq_min freq_max];
yrange = [0 group_max];
domelev_image_handle = imagesc(xrange, yrange, elevation_dom, [0 50]);
set(gca, 'ydir', 'normal', 'Fontsize', 14)
ylabel('Group Range (km)', 'Fontsize', 14)

xlabel('Frequency (MHz)', 'Fontsize',  14)
cb_handle = colorbar('Fontsize', 14);
set(get(cb_handle,'ylabel'), 'string', 'Elevation (degrees)', 'Fontsize', 14, ...
  'verticalalignment','top')


%
% plot the number of hops of the dominant mode 
%
domnumhop_fig_handle = figure('Name', ...
  'Synthetic BSS Ionogram: Dominant mode number of hops', ...
  'NumberTitle', 'off', 'Position', [950 500 750 550]);
xrange = [freq_min freq_max];
yrange = [0 group_max];
domnumhop_image_handle = imagesc(xrange, yrange, num_hops_dom, [0 5]);
set(gca, 'ydir', 'normal', 'Fontsize', 14)

red = [1 0 0];
green = [0 1 0];
blue = [0 0 1];
yellow = [1 1 0];
grey = [0.9 0.9 0.9];
cmap = [grey; blue; green; red; yellow];
colormap(cmap)

ylabel('Group Range (km)', 'Fontsize', 14)
xlabel('Frequency (MHz)', 'Fontsize',  14)
cb_handle = colorbar('Fontsize', 14);
set(get(cb_handle,'ylabel'), 'string', 'Number of Hops', 'Fontsize', 14, ...
  'verticalalignment','top')
set(cb_handle, 'ylim', [1 5], 'ytick', [1.5 2.5 3.5 4.5], ...
    'yticklabel', [1 2 3 4]) 
