clear all
close all
clc

load('/home/shibaji7/PHARLAP_SAMI_INTG/pharlap_4.1.3/sami_parsed.mat');
clear deni denn denn_sum
clear nc ns nz nf nt


alt(end,:) = [];
lat(end,:) = [];
lon(end,:) = [];
dene(end,:,:) = [];

intp.threshold = 65;
intp.modu = 2000;
intp.n = 201;
intp.altv = linspace(0,1000,intp.n);
latr = max(max(lat)) - min(min(lat));
intp.latv = linspace(min(min(lat)) + 1,max(max(lat)) - 1,intp.n);
intp.lon = linspace(min(min(lon)),max(max(lon)),intp.n);
intp.lon = mean(intp.lon);

intp.rangev = get_GC_distance(intp.latv,intp.lon,alt);

intp.alt = meshgrid(intp.altv);
intp.lat = meshgrid(intp.latv);
intp.range = meshgrid(intp.rangev);

%intp.alt(end,:) = [];
%intp.lat(end,:) = [];
%intp.range(end,:) = [];

Ne = dene(:,:,1);
intp.Ne = get_intp_data(alt,lat,Ne,intp);
[intp.Ne,attn] = attnFn(intp.Ne,intp.lat,intp.alt,true);


origin_lat = intp.lat(1,1);
origin_lon = intp.lon;
ray_bear = 180;
irregs_flag = 0;
tol = [1e-8,0.01,10];
nhops = 2;
elevs = 20:5:50;
freqs = 10.*ones(size(elevs));
collision_freq = zeros(size(intp.Ne));
range_inc = intp.range(1,2) - intp.range(1,1);
alt_inc = intp.alt(1,2)-intp.alt(1,1);
height_inc = alt_inc;

irreg = zeros(4,intp.n);
[ray_data, ray_path_data] = ...
    raytrace_2d(origin_lat, origin_lon, elevs, ray_bear, freqs, nhops, ...
    tol, irregs_flag, intp.Ne, intp.Ne, ...
    collision_freq,intp.alt(1,1),alt_inc,range_inc,irreg);

figure
pcolor(intp.latv,intp.altv,intp.Ne)
shading interp
xlabel('latitude')
ylabel('altitude')
colorbar

figure
pfsq_conv = 80.6163849431291e-6;
R12 = 100;
UT = [2001 3 15 7 0];
UT_str = [num2str(UT(3)) '/' num2str(UT(2)) '/' num2str(UT(1)) '  ' ...
    num2str(UT(4), '%2.2d') ':' num2str(UT(5), '%2.2d') 'UT'];
freq_str = [num2str(freqs(1)) 'MHz'];
R12_str = num2str(R12);
lat_str = num2str(origin_lat);
lon_str = num2str(origin_lon);
bearing_str = num2str(ray_bear);
fig_str = [UT_str '   ' freq_str '   R12 = ' R12_str '   lat = ' lat_str ...
    ', lon = ' lon_str ', bearing = ' bearing_str];
set(gcf, 'name', fig_str)
start_range = min(min(intp.range));
end_range = max(max(intp.range));
end_range_idx = fix((end_range-start_range) ./ range_inc) + 1;
start_ht = min(min(intp.alt));
start_ht_idx = 1;
end_ht = max(max(intp.alt));
end_ht_idx = fix(end_ht ./ height_inc) + 1;
iono_pf_subgrid = sqrt(intp.Ne.*pfsq_conv);
%iono_pf_subgrid = intp.Ne;
plot_ray_iono_slice(iono_pf_subgrid, start_range, end_range, range_inc, ...
    start_ht, end_ht, height_inc, ray_path_data, 'color', [1, 1, 0.99], ...
    'linewidth', 2);

%set(gcf,'units','normal')
%pos = get(gcf,'position');
%pos(2) = 0.55;
%set(gcf,'position', pos)
% uncomment the following to print figure to hi-res ecapsulated postscript
% and PNG files
%set(gcf, 'paperorientation', 'portrait')
%set(gcf, 'paperunits', 'cent', 'paperposition', [0 0 61 18])
%set(gcf, 'papertype', 'a4')


sqrt(max(max(intp.Ne))*pfsq_conv)