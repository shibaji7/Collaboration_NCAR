%
% setup general stuff
%
function Eclipse_Raytracing(origin_lat_p,origin_long_p,ray_bear_p, UT_Array_n, dir_path,offset_n)
PHARLAP_OLD_FORMAT_2D='true';
clear t;
%% Default values of the parameters in the dialog box
%UT = [2016 8 21 18 30];        % UT - year, month, day, hour, minute
UT_array = UT_Array_n;
%UT_array = [2016 8 21 17 37;2016 8 21 17 58;2016 8 21 18 14];
%UT_array = [2016 8 21 18 30];
offset = offset_n;
origin_lat = origin_lat_p;           % latitude of the start point of ray
origin_long = origin_long_p;       % longitude of the start point of ray
ray_bear = ray_bear_p;           % Azimuth angle of the rays
end_range = '2000';            % end point of the map
elev_start='15';
elev_inc='5';
elev_end='90'; %45-90,9 deg steps
freq_start='4';               % range of ray frequencies
freq_inc='2';
freq_end='14';

%% The dialog box
% Define default variables
prompt={'origin_lat','origin_long','ray_bear','end_range','freq_start',...
    'freq_inc','freq_end','elev_start','elev_inc','elev_end'};
dlg_title='Input Values';
num_lines=1;
defAns={origin_lat,origin_long,ray_bear,end_range,freq_start,...
    freq_inc,freq_end,elev_start,elev_inc,elev_end};
options='on';
% Take inputs
answer = inputdlg(prompt,dlg_title,num_lines,defAns,options);
% Parse inputs

origin_lat=str2num(answer{1});
origin_long=str2num(answer{2});
ray_bear=str2num(answer{3});
end_range=str2num(answer{4});
freq_start=str2num(answer{5});
freq_inc=str2num(answer{6});
freq_end=str2num(answer{7});
elev_start=str2num(answer{8});
elev_inc=str2num(answer{9});
elev_end=str2num(answer{10});

%% Conversion of the parameters
elevs = elev_start:elev_inc:elev_end;
load('LoadSplinefit.mat');
% UT_end = [2016 8 21 20 10];
R12 = 100;                   % R12 index
speed_of_light = 2.99792458e8;
speed_of_earth = 40070/(24 * 60); %km/hr

%% Other variables declaration
doppler_flag = 1;            % generate ionosphere 5 minutes later so that
% Doppler shift can be calculated
irregs_flag = 0;             % no irregularities - not interested in
% Doppler spread or field aligned irregularities
kp = 0;                      % kp not used as irregs_flag = 0. Set it to a
% dummy value

% generate ionospheric, geomagnetic and irregularity grids
max_range = 5000;      % maximum range for sampling the ionosphere (km)
num_range = 1001;        % number of ranges (must be < 2000)
range_inc = max_range ./ (num_range - 1);  % range cell size (km)
speed_at_lat= cosd(origin_lat) * speed_of_earth ;
range_time_inc=(range_inc/speed_at_lat)/(24*60);

start_height = 0 ;      % start height for ionospheric grid (km)
height_inc = 3;         % height increment (km)
num_heights = 200;      % number of  heights (must be < 2000)

clear iri_options
iri_options.Ne_B0B1_model = 'Bil-2000'; % this is a non-standard setting for
% IRI but is used as an example

tol = 1e-7;          % ODE tolerance
ray = [];
idx = 1;
start_range = 0;
start_range_n = 0;
end_range_n = end_range;
start_range_idx = fix(start_range_n ./ range_inc) + 1;
end_range_idx = fix(end_range_n ./ range_inc) + 1;
mid_idx=1;
start_ht = start_height;
start_ht_idx = 1;
nhops=1; % number of hops. Currently single hop propagation
R12_str = num2str(R12);
lat_str = num2str(origin_lat);
lon_str = num2str(origin_long);
bearing_str = num2str(ray_bear);
end_ht = 597;
end_ht_idx = fix(end_ht ./ height_inc) + 1;
cnt=1;

%% Generate the damping factor matrix
% Each row corresponds to the location of the eclipse from the point of
% origin.
damping_factor=ones(1,num_range);
 for m=1:1
   % mid_km=mid_idx*range_inc;
   mid_idx = fix(offset/range_inc);
%     mid_idx=fix((end_range_idx - start_range_idx)*(m-1)/4 + start_range_idx);
    damping_factor(m,mid_idx)=splineMidTime(0);
    for d=1:mid_idx-1
        damping_factor(m,mid_idx + d)=splineMidTime(d*range_time_inc);
        damping_factor(m,mid_idx - d)=damping_factor(m,mid_idx + d);
    end
    for d=2*mid_idx:num_range
        damping_factor(m,d)=splineMidTime((d-mid_idx)*range_time_inc);
    end
    [x_dam,y_dam] = size(damping_factor(m,:));
    for z=1:y_dam
      if damping_factor(m,z) > 1
      damping_factor(m,z) = 1;
      end
    end
end

%% Invoke PHaRLAP
% for different time of the day, frequencies and elevation angles.
for i=1:size(UT_array,1)
    UT=UT_array(i,:);
    tic
    [orig_iono_pf_grid, iono_pf_grid_5, collision_freq, irreg] = ...
        gen_iono_grid_2d_sneh(origin_lat, origin_long, R12, UT, ray_bear, ...
        max_range, num_range, range_inc, start_height, ...
        height_inc, num_heights, kp, doppler_flag, 'iri2016', ...
        iri_options);
    toc
    for m=1:1
        iono_pf_grid=orig_iono_pf_grid;
        for s=1:num_range % Apply the damping to the ionosphere
            iono_pf_grid(:, s)= iono_pf_grid(:, s) * (damping_factor(m,s))^0.5;
        end
        % convert plasma frequency grid to  electron density in electrons/cm^3
        iono_en_grid = iono_pf_grid.^2 / 80.6164e-6;
        iono_en_grid_5 = iono_pf_grid_5.^2 / 80.6164e-6;
        for freq=freq_start:freq_inc:freq_end
            freqs = freq.*ones(size(elevs));
            %figure(cnt);
            [ray_data, ray_path_data] = ...
                raytrace_2d(origin_lat, origin_long, elevs, ray_bear, freqs, nhops, ...
                tol, irregs_flag, iono_en_grid, iono_en_grid_5, ...
                collision_freq, start_height, height_inc, range_inc, irreg);
            
            UT_str = [num2str(UT(3)) '-' num2str(UT(2)) '-' num2str(UT(1)) '-' ...
                num2str(UT(4), '%2.2d') '-' num2str(UT(5), '%2.2d') 'UT'];
            freq_str = [num2str(freqs(1)) 'MHz'];
            fig_str = [UT_str '   ' freq_str '   R12 = ' R12_str '   lat = ' lat_str ...
                ', lon = ' lon_str ', bearing = ' bearing_str];
            
            set(gcf, 'name', fig_str)
            iono_pf_subgrid = iono_pf_grid(start_ht_idx:end_ht_idx, ...
                start_range_idx:end_range_idx);

            [axis_handle, ray_handle] = plot_ray_iono_slice(iono_pf_subgrid, ...
                start_range, end_range, range_inc, start_ht, end_ht, height_inc, ...
                ray_path_data, 'color', 'w', 'linewidth', 2);

            legend(UT_str, freq_str);
            title(fig_str)
            %figname=strcat(freq_str,'_',UT_str,'_',num2str(m),'.png');
			figname=strcat(freq_str,'-',UT_str,'-',num2str(m));
			t= annotation('textbox',[.43 .52 .4 .5],'String', dir_path);
            t.FontSize = 12;
            t.LineStyle = 'none';
            
            %saveas(gcf,figname)
			set(gcf, 'PaperPositionMode','auto'); 
	        saveas(gcf,[dir_path '\Eclipse_' fig_str '.png']);
            cnt=cnt+1;
        end
    end
end
