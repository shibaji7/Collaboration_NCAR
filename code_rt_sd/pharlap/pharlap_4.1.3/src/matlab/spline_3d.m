%
% Name :
%   spline_3d.m
%
% Purpose :
%   Generates a 3d attenuation function in latitude and longitude.
%   Computation the latitude and and longitude axis to scale 
%
% Calling sequence :
%   spline_3d
%
% Inputs :
% eclipse_centre_lat: latitude for the centre of the eclipse
% eclipse_centre_lon: longitude for the centre of the eclipse
% grid_inc_lat: latitude incerment (deg), same as lat_inc in set_ionogrid_3d for ionospheric grid
% grid_inc_lon: longitude incerment (deg), same as lon_inc in set_ionogrid_3d
%
% Outputs :
% lat_lon_attn: 3d matrix with attenuation facors corresponding to latitude
% and longitude

% input parameters for 3d attenuation function
eclipse_centre_lat = 40;              % latitude for the centre of the eclipse
eclipse_centre_lon = -100;            % longitude for the centre of the eclipse
grid_inc_lat = 1;                     % latitude incerment (deg), same as lat_inc in set_ionogrid_3d for ionospheric grid
grid_inc_lon = 1;                     % longitude incerment (deg), same as lon_inc in set_ionogrid_3d
num_points = 10;                      % to specify number of points on the circle for the cylinder function
rotate_fac = fix(360/num_points) ;    % rotation angle in azimuth for ray2latlon

% 2d attenuation function 
load('LoadSplinefit.mat');      % Load 2d spline fit in time
origin_lat = 38.86;             % Origin Latitude
offset = 5;                     
speed_of_light = 2.99792458e8;
speed_of_earth = 40070/(24 * 60); %km/hr

% generate ionospheric, geomagnetic and irregularity grids
max_range = 5000;       % maximum range for sampling the ionosphere (km)
num_range = 101;        % number of ranges (must be < 2000)
range_inc = max_range ./ (num_range - 1);  % range cell size (km)
speed_at_lat= cosd(origin_lat) * speed_of_earth ;
range_time_inc=(range_inc/speed_at_lat)/(24*60);

start_range = 0;
end_range = 2000;            % end point of the map
start_range_idx = fix(start_range ./ range_inc) + 1;
end_range_idx = fix(end_range ./ range_inc) + 1;
mid_idx=1;

    %% Generate the damping factor matrix for the 3d case from 2d
    % damping_facto corresponds to the location of the eclipse from the point of
    % origin.
    damping_factor=ones(1,num_range);
    m = 1;
    %mid_idx = fix(offset/range_inc);
    damping_factor(m,mid_idx)=splineMidTime(0);
    for d=1:mid_idx-1
        damping_factor(m,mid_idx + d)=splineMidTime(d*range_time_inc);
        damping_factor(m,mid_idx - d)=damping_factor(m,mid_idx + d);
    end
    for d=2*mid_idx:num_range
        damping_factor(m,d)=splineMidTime((d-mid_idx)*range_time_inc);
    end
    [x_dam,y_dam] = size(damping_factor(m,:));
    % attenuation function greater than 1 is set to 1
    for z=1:y_dam
      if damping_factor(m,z) > 1
      damping_factor(m,z) = 1;
      end
    end
    xa = 1:length(damping_factor);
    rad_d = abs((xa-mid_idx)* range_inc);
    %[x_dist,y_dist, z_dist] = cylinder(rad_d);
    % [X,Y,Z] = cylinder(r,n) returns the x-, y-, and z-coordinates of a cylinder based on the profile curve 
    % defined by vector r. The cylinder has n equally spaced points around its circumference.
    [x_dist,y_dist, z_dist] = cylinder(rad_d, num_points);
    x_dist=x_dist+(mid_idx*range_inc);
    y_dist=y_dist+(mid_idx*range_inc);
    % z is the attenuation factor
    [xn,yn] = size(z_dist);
    for l=1:yn
        z_dist(:,l) = damping_factor(m,:);
           
    end	
    figure(8);
    surfc(x_dist,y_dist, z_dist);
    xlabel('x');
    ylabel('y');
    zlabel('attn');
    %% Map the attenuation factor in distance to latitude and longitude
    figure(18)
    % raz2latlon: Converts the ground range and azimuth (from true North) of point from a 
    % particular bore point or origin on the Earth to latitude and longitude.
    for k =1:num_points + 1;
    [latr_mat(:,k),lonr_mat(:,k)] = raz2latlon(rad_d'.*1000, rotate_fac*(k-1), eclipse_centre_lat, eclipse_centre_lon ,'wgs84');
    end
    surfc(latr_mat,lonr_mat, z_dist);
    xlabel('latitude(deg)');
    ylabel('longitude(deg)');
    zlabel('attn');
    
    % store lat, lon and attenuation factor for all points
    for k=1:num_points
    damping_grid(:,:,k)=[latr_mat(:,k),lonr_mat(:,k),z_dist(:,k)];
    end
    
%% Construct a 3 dim matrix of latitude, longitude and corresponding attenuation factor
% concatenate damping_grid for the poinys (num_points) obtained
lat_lon = [damping_grid(:,:,1);damping_grid(:,:,2);damping_grid(:,:,3);damping_grid(:,:,4);damping_grid(:,:,5);damping_grid(:,:,6);damping_grid(:,:,7);damping_grid(:,:,8);damping_grid(:,:,9);damping_grid(:,:,10)]
j=1;
for i=1:length(lat_lon)
    if lat_lon(i,3) ~=1
    lat_lon_mod(j,:) =  lat_lon(i,:);
    j=j+1;
    end
end
[values, order] = sort(lat_lon_mod(:,1));
lat_lon_sort = lat_lon_mod(order,:);

%% Calculation to scale
%run('lat_lon_mod.m');
lat_lon_fix1(:,1)=fix((lat_lon_sort(:,1)*1))/1;
lat_lon_fix1(:,2)=fix((lat_lon_sort(:,2)*1))/1;
lat_lon_fix1(:,3)=lat_lon_sort(:,3);
% remove duplicates
lat_lon_fix_unq1=unique(lat_lon_fix1,'rows');
% find maximum and minimum latitude and longitude values from the matrix
latr_min= min(lat_lon_fix_unq1(:,1));
latr_max= max(lat_lon_fix_unq1(:,1));
lonr_min= min(lat_lon_fix_unq1(:,2));
lonr_max= max(lat_lon_fix_unq1(:,2));
i=1;
j=1;
z_mean = lat_lon_fix_unq1(1,3);
% generate attenuation factors that correspond to the subgrid ionospheric
% grid for the eclipsed region
for latr_num = latr_min:grid_inc_lat:latr_max
for lonr_num = lonr_min:grid_inc_lon:lonr_max
% find required latitude and longitude 
ind = find(lat_lon_fix_unq1(:,1)==latr_num & lat_lon_fix_unq1(:,2)==lonr_num);

if length(ind) ~= 0
% For the same latitude, longitude values, take mean of all the available
% attenuation factor valued
z_mean = mean(lat_lon_fix_unq1(ind,3));
z_std = std(lat_lon_fix_unq1(ind,3));
lat_lon_new(i,:) = [latr_num lonr_num z_mean];
lat_lon_attn(j,:) = [latr_num lonr_num z_mean];
i=i+1;
j=j+1;
else
lat_lon_attn(j,:) = [latr_num lonr_num z_mean];
j=j+1;
end
end
end