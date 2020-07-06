%
% Name :
%   ois_synth_mh.m
%
% Purpose :
%   Synthetic multi-hop OIS ionogram generation using 2d raytracing.
%
% Calling sequence :
%   ois_synth_mh()
%
% Inputs :
%
% Outputs :
%   
% Author:
%   V1.0  M.A. Cervera  13/12/2012
%
%   V1.1  M.A. Cervera  15/02/2013
%      Various Improvements to GUI and bug fixes
%
%   V1.2  M.A. Cervera  06/07/2015
%      Updated the method of calculating the O-X correction - now uses
%      gm_freq_offset.m
%
%   V1.3  M.A. Cervera  02/05/2016
%      Fixed bug (introduced in V1.2) regarding how how the geomagnetic O-X 
%      correction to the 2D raytrace was being applied. Updated to use
%      IRI2016
%
%   V1.4  M.A. Cervera  20/05/2016
%      Updated to use multi-threaded raytrace_2d
%



function [] = ois_synth_mh()

synth_iono = 'Yes';

month_str_init = '1';
hour_str_init = '00';
minute_str_init = '00';
SSN_str_init = '100';
tx_lat_str_init = '-25';
tx_lon_str_init = '135';
rx_lat_str_init = '-10';
rx_lon_str_init = '135';
num_hops_str_init = '4';

while(strcmp(synth_iono, 'Yes'))
  
  %
  % setup general stuff
  %
  prompt = {'Month (1 - 12)', 'UT Hour (00 - 23)', 'UT Minute (00 - 59)', ...
    'Smoothed Sunspot Number (R12 index: 1 - 150)', ...
    'Transmitter Latitude (decimal degrees: -90 to 90)', ...
    'Transmitter Longitude (decimal degrees: -180 to 180)', ...
    'Receiver Latitude (decimal degrees: -90 to 90)', ...
    'Receiver Longitude (decimal degrees: -180 to 180)' ...
    'Number of hops (1 - 10)'};
  dlg_title = 'Input OIS Synth parameters';
  num_lines = 1;
  def = {month_str_init, hour_str_init, minute_str_init, SSN_str_init, ...
         tx_lat_str_init, tx_lon_str_init, rx_lat_str_init, ...
         rx_lon_str_init num_hops_str_init};
  answer = inputdlg(prompt, dlg_title, num_lines, def);
  
  if isempty(answer), return, end
  
  speed_of_light = 2.99792458e8;
  
  month_str  = answer{1};
  hour_str   = answer{2};
  minute_str = answer{3};
  SSN_str    = answer{4};
  tx_lat_str = answer{5};
  tx_lon_str = answer{6};
  rx_lat_str = answer{7};
  rx_lon_str = answer{8};
  num_hops_str = answer{9};
  UT_str = [hour_str ':' minute_str];
  
  UT(1) = 2000;   %year - don't really care what this is as SSN is used
  UT(2) = str2num(month_str);
  UT(3) = str2num(hour_str);
  UT(4) = str2num(minute_str);
  UT(5) = 0;
  switch UT(2)
    case 1, month_str = 'JAN';
    case 2, month_str = 'FEB';
    case 3, month_str = 'MAR';
    case 4, month_str = 'APR';
    case 5, month_str = 'MAY';
    case 6, month_str = 'JUN';
    case 7, month_str = 'JUL';
    case 8, month_str = 'AUG';
    case 9, month_str = 'SEP';
    case 10, month_str = 'OCT';
    case 11, month_str = 'NOV';
    case 12, month_str = 'DEC';
  end
  
  R12 = str2num(SSN_str);
  tx_lat = str2num(tx_lat_str);
  tx_lon = str2num(tx_lon_str);
  rx_lat = str2num(rx_lat_str);
  rx_lon = str2num(rx_lon_str);
  nhops = str2num(num_hops_str);
  
  
  figure(1)
  set(gcf, 'menubar', 'none', 'numbertitle', 'off', 'name', ...
    'Synthetic OIS generation tool:')
  pos = get(gcf,'position');
  pos(4) = 200;
  set(gcf,'position',pos);
  text(0,1,[month_str ' ' UT_str])
  text(0,0.85, ['SSN: ' SSN_str])
  text(0,0.7, ['Tx. lat. ' tx_lat_str])
  text(0,0.6, ['Tx. lon. ' tx_lon_str])
  text(0,0.5, ['Rx. lat. ' rx_lat_str])
  text(0,0.4, ['Rx. lon. ' rx_lon_str])
  set(gca,'visible','off')
  figure(1)
  
  name_str = ['Synthesised OIS :        tx lat ' tx_lat_str ', tx lon ' ...
    tx_lon_str ', rx lat ' rx_lat_str ', rx lon ' rx_lon_str  ...
    '               ' month_str ' ' UT_str '       SSN ' SSN_str];
  
  drawnow
  
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
  max_range = 15000;      % maximum range for sampling the ionosphere (km)
  num_range = 301;        % number of ranges (must be < 2001)
  range_inc = max_range ./ (num_range - 1);  % range cell size (km)
  
  start_height = 60;      % start height for ionospheric grid (km)
  height_inc = 2;         % height increment (km)
  num_heights = 200;      % number of  heights (must be < 2001)
  
  doppler_flag = 0;       % not interested in Doppler shift and spread
  kp = 0;                 % kp not used as doppler_flag = 0. Set it to a 
                          % dummy value 

  h = text(0, 0.2, 'Generating ionosphere...'); drawnow
  figure(1)
  [iono_pf_grid, iono_pf_grid_5, collision_freq, irreg] = ...
       gen_iono_grid_2d(tx_lat, tx_lon, R12, UT, azim_rx, ...
                        max_range, num_range, range_inc, start_height, ...
                        height_inc, num_heights, kp, doppler_flag, 'iri2016');
  set(h, 'string', 'Generating ionosphere... Done.'); drawnow
  
  %convert plasma frequency grid to electron density in electrons/cm^3
  iono_en_grid = iono_pf_grid.^2 / 80.6164e-6;
  iono_en_grid_5 = iono_pf_grid_5.^2 / 80.6164e-6;
  
  %
  % call raytrace - loop over frequency and elevation
  %
  irregs_flag = 0;     % no irregularities
  tol = 1e-7;          % ode solver tolerance  
  
  frequency = [];
  elevation = [];
  group_range= [];
  dev_absorption = [];
  D_Oabsorp = [];
  D_Xabsorp = [];
  fs_loss = [];
  eff_range = [];
  phase_path = [];
  ray_apogee = [];
  ray_apogee_gndr = [];
  plasfrq_at_apogee = [];
  ray_hops = [];
  dgnd_dels = [];
  del_freq_O = [];
  del_freq_X = [];

  colour_str = ['b', 'g', 'r'];
  elev_step = 0.5;
  elevs = [2 : elev_step : 80];
  num_elevs = length(elevs);
  
  
  % Do two passes of raytracing: coarse over the entire frequency range
  % and then fine around the MUF.
  pass_num = 0;
  while(pass_num < 2)
    
    pass_num = pass_num + 1;
    
    if (pass_num == 1)
      start_freq = 2;
      end_freq = 45;
      freq_step = 0.5;
    else
      start_freq = coarse_muf - 0.5;
      end_freq = coarse_muf + 0.5;
      freq_step = 0.05;
    end
    min_freq_stepsize = 0.01;
    freq = start_freq - freq_step;
    old_numrays_to_rx = 0;
    
    % first call to raytrace so pass in the ionospheric and geomagnetic grids
    [ray_data, ray_path_data] = ...
        raytrace_2d(tx_lat, tx_lon, elevs(1), ray_bearing, start_freq, ...
        nhops, tol, irregs_flag, iono_en_grid, iono_en_grid_5, ...
        collision_freq, start_height, height_inc, range_inc, irreg);
    
    % loop over ray frequency
    if (pass_num == 1)
      h1 = text(0, 0.1, 'Synthesizing OIS (first pass):');
    else
      h2 = text(0, 0, 'Synthesizing OIS (second pass):');
    end
    while (freq <= end_freq)
      
      freq = freq + freq_step;
      
      if (pass_num == 1)
        textstr = sprintf('Synthesizing OIS (first pass): freq = %3.1f MHz', ...
		          freq);
        figure(1)
        set(h1,'string',textstr)
      else
        textstr = sprintf('Synthesizing OIS (second pass): freq = %3.1f MHz',...
	                  freq);
        figure(1)
        set(h2,'string',textstr)
      end
      drawnow
      
      % loop over ray elevation 0.5 degree steps
      ray_data = zeros(num_elevs, nhops)*NaN;
      gnd_range = zeros(num_elevs, nhops)*NaN;
      grp_range = zeros(num_elevs, nhops)*NaN;
      labels = zeros(num_elevs, nhops);
      freqs = freq .* ones(size(elevs));
      [ray_data, ray_path_data] = raytrace_2d(tx_lat, tx_lon, elevs, ...
	           ray_bearing, freqs, nhops, tol, irregs_flag);
      for el_idx=1:num_elevs
        for hop_idx = 1:ray_data(el_idx).nhops_attempted
          gnd_range(el_idx, hop_idx) = ray_data(el_idx).ground_range(hop_idx);
          grp_range(el_idx, hop_idx) = ray_data(el_idx).group_range(hop_idx);
          labels(el_idx, hop_idx) = ray_data(el_idx).ray_label(hop_idx);
        end
      end
      
      
      % find the "good rays" i.e. the rays which come to ground OK
      for hop_idx = 1:nhops
        idx_goodray = find(labels(:, hop_idx) == 1);
        if length(idx_goodray) > 3
          
          % Find ray ground ranges which bracket the receiver ground range, do
          % raytracing with finer (0.05 deg) elevation grid within coarse 
          % braketing rays, and finally interpolate to find the ray elevations 
          % and group ranges (and their absorption losses), which will hit
          % the receiver.
          els = elevs(idx_goodray);
          gnd = gnd_range(idx_goodray, hop_idx)';
          grp = grp_range(idx_goodray, hop_idx)';
          dgrp_dels = deriv(grp, els);
          num = length(els);
          grp_to_rx = [];
          
          % loop over all good elevations
          for ii=1:num-1
            
            % find the bracketing rays - ignore those whose rate of change of 
	    % range  with elevation is too large as this indicates we are too
	    % far into a cusp region to be reliable
            if ((gnd(ii) >= range_rx & gnd(ii+1) < range_rx) | ...
                (gnd(ii) <= range_rx & gnd(ii+1) > range_rx)) & ...
                (els(ii+1) - els(ii) < 2*elev_step) & ...
                (abs(dgrp_dels(ii)) < 500) & (abs(dgrp_dels(ii+1)) < 500)
              
              el_step = els(ii+1) - els(ii);
              fine_el_step = el_step ./ 5;
              fine_els = [els(ii) : fine_el_step : els(ii+1)];
              fine_elevs = [];
              fine_gnd = [];
              fine_label = [];
              
              % raytrace at fine elevation steps between bracketing rays
	      freqs = freq .* ones(size(fine_els));
              [ray_data, ray_path_data] = raytrace_2d(tx_lat, tx_lon, ...
		        fine_els, ray_bearing, freqs, hop_idx, tol, irregs_flag);
                
              for idx=1:6
		fine_elev = fine_els(idx);
                if ray_data(idx).nhops_attempted == hop_idx
                  fine_gnd = [fine_gnd ray_data(idx).ground_range(hop_idx)];
                  fine_label = [fine_label ray_data(idx).ray_label(hop_idx)];
                  fine_elevs = [fine_elevs fine_elev];
                end
              end
              
              % interpolate to get elevation to launch ray to hit rx and
              % raytrace at this elevation to get all the other required
              % quantities 
              if (isempty(find(fine_label < 1)) & length(fine_label >=3))
                elev_torx = interp1(fine_gnd, fine_elevs, range_rx, 'pchip');
                
                [ray_data, ray_path_data] = raytrace_2d(tx_lat, tx_lon, ...
		      elev_torx, ray_bearing, freq, hop_idx, tol, irregs_flag);
                
                if ray_data.ray_label == 1
                  elevation = [elevation elev_torx];
                  group_range = [group_range ray_data.group_range(hop_idx)];
                  phase_path =  [phase_path ray_data.phase_path(hop_idx)];
                  dev_absorption = ...
		      [dev_absorption ray_data.deviative_absorption(hop_idx)];
                  eff_range = ...
		      [eff_range ray_data.effective_range(hop_idx)];
                
                  gnd_fs_loss = 0;
                  O_absorp = 0;
                  X_absorp = 0;
                  for kk = 1:hop_idx
                    ray_apogee  = ray_data.apogee(kk);
                    ray_apogee_gndr  = ray_data.gnd_rng_to_apogee(kk);
                    [ray_apogee_lat, ray_apogee_lon] = raz2latlon( ...
                        ray_apogee_gndr, ray_bearing, tx_lat, tx_lon, 'wgs84');
		    plasfrq_at_apogee = ray_data.plasma_freq_at_apogee(kk);
                    if kk == 1
		      % calculate geo-mag splitting factor - assume that it
		      % is the same for all hops (really need to calculate
		      % separately for each hop)
  	              [del_fo, del_fx] = ...
	                    gm_freq_offset(ray_apogee_lat, ray_apogee_lon, ...
                                           ray_apogee, ray_bearing, ...
				           freq, plasfrq_at_apogee, UT);
                    end
                    O_absorp = O_absorp + abso_bg(ray_apogee_lat, ...
			ray_apogee_lon, elev_torx, freq + del_fo, ...
			UT, R12, 1);
                    X_absorp = X_absorp + abso_bg(ray_apogee_lat, ...
			ray_apogee_lon, elev_torx, freq + del_fx, ...
			UT, R12, 0);

                    if kk > 1
                      fs_lat = ray_data.lat(kk-1);
                      fs_lon = ray_data.lon(kk-1);
                      gnd_fs_loss = gnd_fs_loss + ...
                        ground_fs_loss(fs_lat, fs_lon, elev_torx, freq);
                    end
                  end
                
                  D_Oabsorp = [D_Oabsorp O_absorp];
                  D_Xabsorp = [D_Xabsorp X_absorp];
                  fs_loss = [fs_loss gnd_fs_loss];
	          del_freq_O = [del_freq_O, del_fo];
	          del_freq_X = [del_freq_X, del_fx];
		
                  frequency = [frequency freq];
                  ray_hops = [ray_hops hop_idx];
		end  % of if label == 1
		
              end
              
            end   % of find bracketing rays
            
          end   % of loop over "good" elevations
          
        end   % of "if length(idx_goodray) > 3"
      end   % of loop over nhops
      
    end   % of frequency loop
    
    coarse_muf = max(frequency);
    
    % don't do second pass if coarse_muf is not a real number (ie if there
    % is no propagtion)
    no_propagation = 0;
    if isempty(coarse_muf) | isnan(coarse_muf) | isinf(coarse_muf) | ...
        ~isreal(coarse_muf) | ~isnumeric(coarse_muf)
      pass_num = 3;
      no_propagation = 1;
    end
    
  end  % of pass_num
  
  delete(gcf)
  
  
  %
  % If there is a propagation mode then plot the OIS
  %
  if no_propagation
    
    mh = msgbox(['No propagation modes available for the number of hops ' ...
            'specified. Try again with more hops.'])
    waitfor(mh)
    
    month_str_init = num2str(UT(2));
    hour_str_init = hour_str;
    minute_str_init = minute_str;
    SSN_str_init = SSN_str;
    tx_lat_str_init = tx_lat_str;
    tx_lon_str_init = tx_lon_str;
    rx_lat_str_init = rx_lat_str;
    rx_lon_str_init = rx_lon_str;
    num_hops_str_init = num_hops_str;
    
  else
    
    %
    % Apply O-X correction
    %
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
    power_O = 10*log10(pow) - dev_absorption - D_Oabsorp - fs_loss;
    power_X = 10*log10(pow) - dev_absorption - D_Xabsorp - fs_loss;


    %
    % plot the synthesized oblique ionogram
    %

    % first power
    figure('Name', 'Close this window to continue', 'NumberTitle', 'off', ...
      'Position', [200 200 1200 450]);
    subplot(1,3,1)
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

      h = plot(frequency(ii), group_range(ii), 'o', 'markerfacecolor', ...
	        m_colour_O, 'markeredgecolor', m_colour_O, 'markersize', 5);
      set(h, 'zdata', power_O(ii)); 
      hold on

    end
    hold off

    min_grp = fix(min(group_range_O)/100) * 100;
    max_grp = round(max(group_range_X)/100 + 0.5) * 100;
    max_grp = max([max_grp min_grp+500]);
    set(gca,'ylim',[min_grp max_grp])
    set(gca,'fontsize',12)
    ylabel('Group Range (km)', 'fontsize', 12)
    xlabel('Frequency (MHz)', 'fontsize', 12)
    c_tick_step = 5;
    c_tickval = [0 : c_tick_step/(max_barrange-min_barrange) : 1];
    if verLessThan('matlab', '8.4.0') c_tickval = c_tickval*num_colours; end
    c_ticklabel = num2str([min_barrange:c_tick_step:max_barrange]');
    cb = colorbar('ytick', c_tickval, 'yticklabel', c_ticklabel);
    set(get(cb, 'ylabel'), 'String', 'Power (dBW)', 'fontsize', 12)
    set(gca,'position', [0.07 0.12 0.2 0.81])


    % plot the elevation
    subplot(1,3,2)
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

      h = plot(frequency(ii), group_range(ii), 'o', 'markerfacecolor', ...
	       m_colour_O, 'markeredgecolor', m_colour_O, 'markersize', 5);
      set(h, 'zdata', elevation(ii));
      hold on

    end
    hold off

    set(gca,'ylim',[min_grp max_grp])
    set(gca,'fontsize',12)
    title([month_str ' ' UT_str '       SSN ' SSN_str ...
      '     1 Watt radiator, iostropic antennas       tx lat ' ...
      tx_lat_str ', tx lon ' tx_lon_str ', rx lat ' rx_lat_str ...
      ', rx lon ' rx_lon_str], 'fontsize', 12);
    ylabel('Group Range (km)', 'fontsize', 12)
    xlabel('Frequency (MHz)', 'fontsize', 12)
    c_tick_step = 10;
    c_tickval = [0 : c_tick_step/(max_barrange-min_barrange) : 1];
    if verLessThan('matlab', '8.4.0') c_tickval = c_tickval*num_colours; end
    c_ticklabel = num2str([min_barrange:c_tick_step:max_barrange]');
    cb = colorbar('ytick', c_tickval, 'yticklabel', c_ticklabel);
    set(get(cb, 'ylabel'), 'String', 'Elevation (degrees)', 'fontsize', 12)
    set(gca,'position', [0.425 0.12 0.2 0.81])


    % plot the number of hops
    subplot(1,3,3)
    hold off
    hop_handles = zeros(1,nhops)*NaN;
    for ii = 1:length(frequency)
      switch ray_hops(ii)
        case 1
          colour_str = 'b';
          edge_colour_str = 'b';
        case 2
          colour_str = 'g';
          edge_colour_str = 'g';
        case 3
          colour_str = 'r';
          edge_colour_str = 'r';
        case 4
          colour_str = 'y';
          edge_colour_str = 'k';
        case 5
          colour_str = 'k';
          edge_colour_str = 'k';
        case 6
          colour_str = 'm';
          edge_colour_str = 'm';
        case 7
          colour_str = 'c';
          edge_colour_str = 'c';
        case 8
          colour_str = [0.95 0.6 0];
          edge_colour_str = [0.95 0.6 0];
        case 9
          colour_str = [0.5 0.5 0.5];
          edge_colour_str = [0.5 0.5 0.5];
        case 10
          colour_str = [0.65 0.65 0.65];
          edge_colour_str = [0.65 0.65 0.65];
      end

      h = plot(frequency(ii), group_range(ii), 'o', 'markerfacecolor', ...
        colour_str, 'markeredgecolor', edge_colour_str, 'markersize', 5);
      hold on

      hop_handles(ray_hops(ii)) = h;
    end
    hold off

    set(gca,'ylim',[min_grp max_grp])
    set(gca,'fontsize',12)

    ylabel('Group Range (km)', 'fontsize', 12)
    xlabel('Frequency (MHz)', 'fontsize', 12)

    legend_labels = [];
    for ii = 1:nhops
      if ~isnan(hop_handles(ii))
        legend_labels = [legend_labels {[num2str(ii) '-hop']}];
      end
    end

    idx_notnan = find(~isnan(hop_handles));
    legend(hop_handles(idx_notnan), legend_labels, 'fontsize', 12)

    set(gca,'position', [0.77 0.12 0.2 0.81])


    set(gcf, 'paperunits', 'cent', 'papertype', 'A4', 'paperorientation', ...
	'landscape', 'PaperPosition', [0 8 29.5 10])

    waitfor(gcf)

    % synthesize another OIS ?
    res = questdlg('Do you want to synthesize another OIS?', 'title', ...
	           'Yes', 'No', 'Yes');
    synth_iono = res;

    delete(gcf)

    month_str_init = num2str(UT(2));
    hour_str_init = hour_str;
    minute_str_init = minute_str;
    SSN_str_init = SSN_str;
    tx_lat_str_init = tx_lat_str;
    tx_lon_str_init = tx_lon_str;
    rx_lat_str_init = rx_lat_str;
    rx_lon_str_init = rx_lon_str;
    num_hops_str_init = num_hops_str;

  end    % of "if no_propagation"
  
end   % of "while(strcmp(synth_iono, 'Yes'))"


end  % of function
