% script for comparing raytrace plots of normal and eclipsed ionosphere for fixed elevation,
% varing azimuths
close all
cnt = 1;

for freq=freq_start:freq_inc:freq_end
    freqs = freq.*ones(size(ray_bears)); % frequency (MHz)
fprintf('Generating %d ''no-field'' rays ...', num_raybears);
%
% call raytrace
%
tic
[ray_data_N, ray_N, ray_state_vec_N] = ...
  raytrace_3d(origin_lat, origin_long, origin_ht, elevs, ray_bears, freqs, ...
              OX_mode, nhops, tol, orig_iono_en_grid, orig_iono_en_grid_5, ...
	          collision_freq, iono_grid_parms, Bx, By, Bz, ...
	          geomag_grid_parms);
          
  [Eclipse_ray_data_N, Eclipse_ray_N, Eclipse_ray_state_vec_N] = ...
  raytrace_3d(origin_lat, origin_long, origin_ht, elevs, ray_bears, freqs, ...
              OX_mode, nhops, tol, iono_en_grid, iono_en_grid_5, ...
	          collision_freq, iono_grid_parms, Bx, By, Bz, ...
	          geomag_grid_parms);

toc
UT_str = [num2str(UT(3)) '-' num2str(UT(2)) '-' num2str(UT(1)) ', ' ...
            num2str(UT(4), '%2.2d') ':' num2str(UT(5), '%2.2d') ' UT'];
        freq_str = [num2str(freqs(1)) 'MHz'];
        fig_str = [UT_str ' ','  ' freq_str];
        
        %set(gcf, 'name', fig_str)
figure(cnt)
for rayId=1:num_raybears
  num = length(ray_N(rayId).lat);
  ground_range = zeros(1, num);
  lat = ray_N(rayId).lat;
  lon = ray_N(rayId).lon;    
  ground_range(2:num) = latlon2raz(lat(2:num), lon(2:num), origin_lat, ...
      origin_long,'wgs84')/1000.0;
  ray_N(rayId).ground_range = ground_range;
end

fprintf('\n')

% plot the rays

for ii = 1:1:num_raybears
  plot3(Eclipse_ray_N(ii).lat,  Eclipse_ray_N(ii).lon, Eclipse_ray_N(ii).height, 'r', 'LineWidth', 1);
  hold on
  plot3(ray_N(ii).lat,  ray_N(ii).lon, ray_N(ii).height, '--b', 'LineWidth', 1)
  hold on
end  
legend('Eclipse' , 'No Eclipse','Location','northeast');
legend('boxoff');
hold off
grid on
title([UT_str ', ' freq_str ', elevation:' num2str(elevs(1))]);
xlabel('latitude (deg)')
ylabel('longitude (deg)')
zlabel('Height (km)')

% view(90,-90) top view
% view(90,0) lon vs height
cnt=cnt+1;
end