%
% Name :
%   ois_synth.m
%
% Purpose :
%   Synthetic single hop OIS ionogram generation using 2d raytracing.
%
% Calling sequence :
%   ois_synth
%
% Inputs :
%   
%
% Outputs :
%   
%
% Author:
%   V1.0  M.A. Cervera  19/06/2006
%
%   V1.1  M.A. Cervera  28/09/2006  
%      Minor update against raytrace_2d
%
%   V2.0  M.A. Cervera  15/06/2007  
%      Updated for Pharlap V2 and includes the effect of different D-region 
%      absorptions on O and X mode rays calulated from the new D region 
%      absorption routine.
%
%   V2.1 M.A. Cervera  12/03/2008
%      Minor update against raytrace_2d. Minor improvements in the alogrithm
%      which finds the ground range bracketing rays.
%
%   V2.2 M.A. Cervera  01/05/2008
%      Modified to work with updated raytrace_2d (for pre-release version of 
%      PHaRLAP 3.0)
%   
%   V2.3 M.A. Cervera  11/09/2009
%      Corrected error with calculation of D region absorption
%
%   V2.4  M.A. Cervera  19/05/2011
%      More efficient handling of ionospheric  and geomagnetic grids grids in
%      call to raytrace_3d  
%
%   V2.5  M.A. Cervera  04/05/2012
%      Fixed bug in power calculation in the radar equation (incorrect 
%      radio wave frequency was used).
%
%   V2.6  M.A. Cervera  06/07/2015
%      Updated the method of calculating the O-X correction - now uses
%      gm_freq_offset.m
%
%   V2.7  M.A. Cervera  02/05/2016
%      Fixed bug (introduced in V2.7) regarding how how the geomagnetic O-X 
%      correction to the 2D raytrace was being applied. Updated to use
%      IRI2016
%
%   V2.6  M.A. Cervera  20/05/2016
%      Updated to use multi-threaded raytrace_2d
%

%
% setup general stuff
%
UT = [2000 1 15 8 0];         % UT - year, month, day, hour, minute
R12 = 116;                    % yearly smoothed sunspot number
speed_of_light = 2.99792458e8;

tx_lat = -12.62;              % latitude of the start point of ray (Darwin)
tx_lon = 131.25;              % longitude of the start point of ray (Darwin)
rx_lat = -23.5230;            % latitude of the end point of ray
rx_lon = 133.678;             % longitude of the end point of ray


% obtain ground range and azimuth of receiver from transmitter
[range_rx, azim_rx] = latlon2raz(rx_lat, rx_lon, tx_lat, tx_lon);
ray_bearing = azim_rx;        % assume omni-directional antnenna => no coning 
range_rx = range_rx / 1000.0; % range now in km

tx_pow = 1;                   % 1 Watt transmitter power
gain_tx = 1.0;                % isotropic transmit antenna
gain_rx = 1.0;                % isotropic receive antenna

%
% generate ionospheric and irregularity grids
%
max_range = 2000;      % maximum range for sampling the ionosphere (km)
num_range = 201;        % number of ranges (must be < 2001)
range_inc = max_range ./ (num_range - 1);  % range cell size (km)

start_height = 60;      % start height for ionospheric grid (km)
height_inc = 2;         % height increment (km)
num_heights = 200;      % number of  heights (must be < 2001)

doppler_flag = 0;            % not interested in Doppler shift and spread
kp = 0;                      % kp not used as doppler_flag = 0. Set it to a 
                             % dummy value 
			     
tic
[iono_pf_grid, iono_pf_grid_5, collision_freq, irreg] = ...
    gen_iono_grid_2d(tx_lat, tx_lon, R12, UT, azim_rx, ...
                     max_range, num_range, range_inc, start_height, ...
		     height_inc, num_heights, kp, doppler_flag);
toc
%convert plasma frequency grid to electron density in electrons/cm^3
iono_en_grid = iono_pf_grid.^2 / 80.6164e-6;
iono_en_grid_5 = iono_pf_grid_5.^2 / 80.6164e-6;

%
% call raytrace - loop over frequency and elevation
%
irregs_flag = 0;     % no irregularities
nhops = 1;           % number of hops
tol = 1e-7;          % ode solver tolerance


% loop over ray frequency, starting with 0.5 MHz steps - decrease step size 
% aroung the nose
frequency = [];
elevation = [];
group_range= [];
dev_absorption = [];
D_Oabsorp = [];
D_Xabsorp = [];
eff_range = [];
phase_path = [];
ray_apogee = [];
ray_apogee_gndr = [];
dgnd_dels = [];
ox_correction = [];
del_freq_O = [];
del_freq_X = [];

elev_step = 0.5;
elevs = [2 : elev_step : 80];
num_elevs = length(elevs);
start_freq = 2;
end_freq = 100;
freq_step = 0.5;
min_freq_stepsize = 0.01;
freq = start_freq - freq_step;
old_numrays_to_rx = 0;
tic

% first call to raytrace so pass in the ionospheric and geomagnetic grids 
[ray_data, ray_path_data] = ...
    raytrace_2d(tx_lat, tx_lon, elevs(1), ray_bearing, start_freq, ...
	       nhops, tol, irregs_flag, iono_en_grid, iono_en_grid_5, ...
	       collision_freq, start_height, height_inc, range_inc, irreg);

while (freq <= end_freq) 
 
  freq = freq + freq_step;

  % call raytrace_2D for a fan of rays at a single frequency: ray elevations at
  % 0.5 degree steps 
  elevs = [2 : elev_step : 80];    % reset the elevs array for current freq
  num_elevs = length(elevs);
  ray_data = zeros(1,num_elevs).*NaN;
  gnd_range = zeros(1,num_elevs).*NaN;
  grp_range = zeros(1,num_elevs).*NaN;
  labels = zeros(1,num_elevs).*NaN;  
  freqs = freq .* ones(size(elevs));
  [ray_data, ray_path_data] = ...
   	    raytrace_2d(tx_lat, tx_lon, elevs, ray_bearing, freqs, nhops, ...
	                tol, irregs_flag);
 
  for idx = 1:length(elevs)
    gnd_range(idx) = ray_data(idx).ground_range;
    grp_range(idx) = ray_data(idx).group_range;
    labels(idx) = ray_data(idx).ray_label;
  end

  
  % find the "good rays" i.e. the rays which come to ground OK
  idx_goodray = find(labels == 1);
  if length(idx_goodray) > 3

    % Find ray ground ranges which bracket the receiver ground range, do
    % raytracing with finer (0.05 deg) elevation grid within coarse braketing
    % rays, and finally interpolate to find the ray elevations and group ranges
    % (and their absorption losses), which will hit the receiver. 
    els = elevs(idx_goodray);
    gnd = gnd_range(idx_goodray);
    grp = grp_range(idx_goodray);
    dgrp_dels = deriv(grp, els);
    num = length(els);
    grp_to_rx = [];

    % loop over all good elevations
    for ii=1:num-1
      
      % find the bracketing rays - ignore those whose rate of change of range
      % with elevation is too large as this indicates we are too far into a cusp
      % region to be reliable
      if ((gnd(ii) >= range_rx & gnd(ii+1) < range_rx) | ...
	  (gnd(ii) <= range_rx & gnd(ii+1) > range_rx)) & ...
          (els(ii+1) - els(ii) < 2*elev_step) & ...
          (abs(dgrp_dels(ii)) < 500) & (abs(dgrp_dels(ii+1)) < 500)

	el_step = els(ii+1) - els(ii);
	fine_el_step = el_step ./ 5;
	fine_elevs = [els(ii) : fine_el_step : els(ii+1)];
	fine_gnd = [];
	fine_grp = [];
	fine_phase = [];
	fine_dev_absorp = [];    % this is only deviative absorption
	fine_eff_range = [];
	fine_ray_apogee = [];
	fine_ray_apogee_gndr = [];
	fine_plasfrq_at_apogee = [];
	fine_label = [];
	freqs = freq .* ones(size(fine_elevs));
	[ray_data, ray_path_data] = ...
	         raytrace_2d(tx_lat, tx_lon, fine_elevs, ray_bearing, ...
		             freqs, nhops, tol, irregs_flag);
        for fine_el_idx=1:6     
	  fine_gnd = [fine_gnd ray_data(fine_el_idx).ground_range];
	  fine_grp = [fine_grp ray_data(fine_el_idx).group_range];
	  fine_phase = [fine_phase ray_data(fine_el_idx).phase_path];
	  fine_dev_absorp = ...
	      [fine_dev_absorp ray_data(fine_el_idx).deviative_absorption];
	  fine_eff_range = ...
	      [fine_eff_range ray_data(fine_el_idx).effective_range];
	  fine_ray_apogee = [fine_ray_apogee ray_data(fine_el_idx).apogee];
	  fine_ray_apogee_gndr = ...
	      [fine_ray_apogee_gndr ray_data(fine_el_idx).gnd_rng_to_apogee];
	  fine_plasfrq_at_apogee = [fine_plasfrq_at_apogee ...
		ray_data(fine_el_idx).plasma_freq_at_apogee];
	  fine_label = [fine_label, ray_data(fine_el_idx).ray_label];
       
	end

	if isempty(find(fine_label < 1))
	  fine_dgnd_dels = deriv(fine_gnd, fine_elevs);

	  elev_torx = interp1(fine_gnd, fine_elevs, range_rx, 'pchip');
	  elevation = [elevation elev_torx];   
	  grp_torx = interp1(fine_elevs, fine_grp, elev_torx, 'pchip');
	  grp_to_rx = [grp_to_rx grp_torx];
	  group_range = [group_range grp_torx];
	  phase_torx = interp1(fine_elevs, fine_phase, elev_torx, 'pchip');
	  phase_path =  [phase_path phase_torx];
	  dev_abso_torx = interp1(fine_elevs, fine_dev_absorp, elev_torx, ...
	      'pchip');
	  dev_absorption = [dev_absorption dev_abso_torx];
	  eff_range_torx = interp1(fine_elevs, fine_eff_range, elev_torx, ...
	      'pchip');
	  eff_range = [eff_range eff_range_torx];
	  ray_apogee_torx = interp1(fine_elevs, fine_ray_apogee, elev_torx, ...
	      'pchip'); 
	  ray_apogee  = [ray_apogee, ray_apogee_torx];
	  ray_apogee_gndr_torx = interp1(fine_elevs, fine_ray_apogee_gndr, ...
	      elev_torx, 'pchip'); 
	  ray_apogee_gndr  = [ray_apogee_gndr, ray_apogee_gndr_torx];
	  dgnd_dels_torx = interp1(fine_elevs, fine_dgnd_dels, elev_torx, ...
	      'pchip'); 
	  dgnd_dels = [dgnd_dels, dgnd_dels_torx];

 	  [ray_apogee_lat, ray_apogee_lon] = raz2latlon( ...
 	        ray_apogee_gndr_torx, ray_bearing, tx_lat, tx_lon, 'wgs84');
	  
          plas_freq_at_apogee = ...
	          interp1(fine_elevs, fine_plasfrq_at_apogee, ...
		          elev_torx, 'pchip');
	  [del_fo_torx, del_fx_torx] = ...
	          gm_freq_offset(ray_apogee_lat, ray_apogee_lon, ...
                                 ray_apogee_torx, ray_bearing, freq, ...
				 plas_freq_at_apogee, UT);
	  del_freq_O = [del_freq_O, del_fo_torx];
	  del_freq_X = [del_freq_X, del_fx_torx];
			     
	  O_absorp = abso_bg(ray_apogee_lat, ray_apogee_lon, elev_torx, ...
	 	freq + del_fo_torx, UT, R12, 1);
	  X_absorp = abso_bg(ray_apogee_lat, ray_apogee_lon, elev_torx, ...
		freq + del_fx_torx, UT, R12, 0);
	  D_Oabsorp = [D_Oabsorp O_absorp];
	  D_Xabsorp = [D_Xabsorp X_absorp];

	  frequency = [frequency ones(size(elev_torx))*freq];

	end

      end

    end


    % If we are near a nose then decrease the frequency step size. Once a
    % minimum step size has been reached then continue with the usual stepsize
    % past the nose.  
    numrays_to_rx = length(grp_to_rx);
    if old_numrays_to_rx > numrays_to_rx
      if ((freq_step ./ 2) > min_freq_stepsize)
	freq = freq - freq_step;
	freq_step = freq_step ./ 5;
      else
	freq_step = 0.5;
	old_numrays_to_rx = 0; 

	% If we are past the F2 nose there is no more propagation so set freq
	% to be greater then end_freq. This will quit the synthesis.
	if (numrays_to_rx == 0) freq = end_freq + 1; end
      end
    else
      old_numrays_to_rx = numrays_to_rx; 
    end


    % plot this group_range solution on the group vs frequency plot
    if numrays_to_rx > 0
      figure(3)
      plot(ones(size(grp_to_rx))*freq, grp_to_rx, 'b.')
      hold on
      mingrp = fix(range_rx ./ 100) .* 100;
      maxgrp = mingrp + 1000;
      set(gca, 'ylim', [mingrp maxgrp],'xlim',[1,35])
      drawnow
    end  
  
  end   % of "if length(idx_goodray) > 3"

end   % of frequency loop
toc
figure(3)
hold off


%
% Apply O-X correction
%
% group_correction = abs(group_range - phase_path) .* ox_correction ./ freqs;
% freqs_O = freqs - ox_correction;
% freqs_X = freqs + ox_correction;
% group_range_O = group_range - group_correction;
% group_range_X = group_range + group_correction;
frequency_O = frequency + del_freq_O;
frequency_X = frequency + del_freq_X;
del_group_O = (group_range - phase_path) .* del_freq_O ./ frequency;
del_group_X = (group_range - phase_path) .* del_freq_X ./ frequency;
group_range_O = group_range + del_group_O;
group_range_X = group_range + del_group_X;


%
% power calculations
%

% one-way radar equation
wavelen = speed_of_light ./ (frequency .* 1e6);
pow = tx_pow * gain_tx .* gain_rx .* (wavelen.^2 ./ (4.*pi)) ./ ...
         (4.*pi .* eff_range.^2);

% ionospheric absorption terms for O and X modes
power_O = 10*log10(pow) - dev_absorption - D_Oabsorp;
power_X = 10*log10(pow) - dev_absorption - D_Xabsorp;


%
% plot the synthesized oblique ionogram
%

% first power
figure(4)
cmap = colormap;
num_colours = length(cmap);

max_barrange = -110;
min_barrange = -160;
scale_factor = (num_colours - 1) ./ (max_barrange - min_barrange);

hold off
for ii = 1:length(frequency)
  colour_idx = round(scale_factor .* (power_O(ii) - min_barrange)) + 1;
  if colour_idx < 1 colour_idx = 1; end
  if colour_idx > length(cmap) colour_idx = length(cmap); end
  m_colour_O = cmap(colour_idx, :);
  
  colour_idx = round(scale_factor .* (power_X(ii) - min_barrange)) + 1;
  if colour_idx < 1 colour_idx = 1; end
  if colour_idx > length(cmap) colour_idx = length(cmap); end
  m_colour_X = cmap(colour_idx, :);

  plot(frequency_O(ii), group_range_O(ii), 'o', 'markerfacecolor', ...
       m_colour_O, 'markeredgecolor', m_colour_O, 'markersize', 5);
  plot(frequency_X(ii), group_range_X(ii), 'o', 'markerfacecolor', ...
       m_colour_X, 'markeredgecolor', m_colour_X, 'markersize', 5);
  hold on
  
end
hold off

min_grp = fix(min(group_range_O)/100) * 100;
max_grp = round(max(group_range_X)/100 + 0.5) * 100;
max_grp = max([max_grp min_grp+500]);
set(gca,'ylim',[min_grp max_grp])
set(gca,'fontsize',18)

title('1 Watt radiator, iostropic antennas')
ylabel('Group Range (km)', 'fontsize', 18)
xlabel('Frequency (MHz)', 'fontsize', 18)
c_tick_step = 5;
c_tickval = [0 : c_tick_step/(max_barrange-min_barrange) : 1];
if verLessThan('matlab', '8.4.0') c_tickval = c_tickval*num_colours; end
c_ticklabel = num2str([min_barrange:c_tick_step:max_barrange]');
cb = colorbar('ytick', c_tickval, 'yticklabel', c_ticklabel);
set(get(cb, 'ylabel'), 'String', 'Power (dBW)', 'fontsize', 18)


% plot the elevation
figure(5)
max_barrange = 60;
min_barrange = 0;
scale_factor = (num_colours - 1) ./ (max_barrange - min_barrange);

hold off
for ii = 1:length(frequency)
  colour_idx = round(scale_factor .* (elevation(ii) - min_barrange)) + 1;
  if colour_idx < 1 colour_idx = 1; end
  if colour_idx > length(cmap) colour_idx = length(cmap); end
  m_colour_O = cmap(colour_idx, :);
  
  colour_idx = round(scale_factor .* (elevation(ii) - min_barrange)) + 1;
  if colour_idx < 1 colour_idx = 1; end
  if colour_idx > length(cmap) colour_idx = length(cmap); end
  m_colour_X = cmap(colour_idx, :);

  plot(frequency_O(ii), group_range_O(ii), 'o', 'markerfacecolor', ...
       m_colour_O, 'markeredgecolor', m_colour_O, 'markersize', 5);
  plot(frequency_X(ii), group_range_X(ii), 'o', 'markerfacecolor', ...
       m_colour_X, 'markeredgecolor', m_colour_X, 'markersize', 5);
  hold on
  
end
hold off

set(gca,'ylim',[min_grp max_grp])
set(gca,'fontsize',18)

ylabel('Group Range (km)', 'fontsize', 18)
xlabel('Frequency (MHz)', 'fontsize', 18)
c_tick_step = 10;
c_tickval = [0 : c_tick_step/(max_barrange-min_barrange) : 1];
if verLessThan('matlab', '8.4.0') c_tickval = c_tickval*num_colours; end
c_ticklabel = num2str([min_barrange:c_tick_step:max_barrange]');
cb = colorbar('ytick', c_tickval, 'yticklabel', c_ticklabel);
set(get(cb, 'ylabel'), 'String', 'Elevation (degrees)', 'fontsize', 18)


% plot the ray apogee
figure(6)
max_barrange = 500;
min_barrange = 0;
scale_factor = (num_colours - 1) ./ (max_barrange - min_barrange);

hold off
for ii = 1:length(frequency)
  colour_idx = round(scale_factor .* (ray_apogee(ii) - min_barrange)) + 1;
  if colour_idx < 1 colour_idx = 1; end
  if colour_idx > length(cmap) colour_idx = length(cmap); end
  m_colour_O = cmap(colour_idx, :);
  
  colour_idx = round(scale_factor .* (ray_apogee(ii) - min_barrange)) + 1;
  if colour_idx < 1 colour_idx = 1; end
  if colour_idx > length(cmap) colour_idx = length(cmap); end
  m_colour_X = cmap(colour_idx, :);

  plot(frequency_O(ii), group_range_O(ii), 'o', 'markerfacecolor', ...
       m_colour_O, 'markeredgecolor', m_colour_O, 'markersize', 5);
  plot(frequency_X(ii), group_range_X(ii), 'o', 'markerfacecolor', ...
       m_colour_X, 'markeredgecolor', m_colour_X, 'markersize', 5);
  hold on
  
end
hold off

set(gca,'ylim',[min_grp max_grp])
set(gca,'fontsize',18)

ylabel('Group Range (km)', 'fontsize', 18)
xlabel('Frequency (MHz)', 'fontsize', 18)
c_tick_step = 50;
c_tickval = [0 : c_tick_step/(max_barrange-min_barrange) : 1];
if verLessThan('matlab', '8.4.0') c_tickval = c_tickval*num_colours; end
c_ticklabel = num2str([min_barrange:c_tick_step:max_barrange]');
cb = colorbar('ytick', c_tickval, 'yticklabel', c_ticklabel);
set(get(cb, 'ylabel'), 'String', 'Ray Apogee (km)', 'fontsize', 18)
