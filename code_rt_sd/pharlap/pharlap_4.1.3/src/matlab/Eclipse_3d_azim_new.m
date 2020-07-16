%
% Name :
%   Eclipse_3d_azim_new.m adapted from ray_test_3d.m
%
% Purpose :
%  script for raytracing for eclipsed ionosphere
%
% Calling sequence :
%   Eclipse_3d_azim_new
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

for freq=freq_start:freq_inc:freq_end
    freqs = freq.*ones(size(ray_bears));
fprintf('Generating %d ''no-field'' rays ...', num_raybears);
tic
[Eclipse_ray_data_N, Eclipse_ray_N, Eclipse_ray_state_vec_N] = ...
  raytrace_3d(origin_lat, origin_long, origin_ht, elevs, ray_bears, freqs, ...
              OX_mode, nhops, tol, iono_en_grid, iono_en_grid_5, ...
	          collision_freq, iono_grid_parms, Bx, By, Bz, ...
	          geomag_grid_parms);
toc
UT_str_eclipse = [num2str(UT(3)) '-' num2str(UT(2)) '-' num2str(UT(1)) '-' ...
            num2str(UT(4), '%2.2d') '-' num2str(UT(5), '%2.2d') 'UT'];
        freq_str_eclipse = [num2str(freqs(1)) 'MHz' 'Eclipse'];
        fig_str_eclipse = [UT_str_eclipse '   ' freq_str_eclipse];
        
        %set(gcf, 'name', fig_str)
figure(cnt)
for rayId=1:num_raybears
  num = length(Eclipse_ray_N(rayId).lat);
  ground_range = zeros(1, num);
  lat = Eclipse_ray_N(rayId).lat;
  lon = Eclipse_ray_N(rayId).lon;    
  ground_range(2:num) = latlon2raz(lat(2:num), lon(2:num), origin_lat, ...
      origin_long,'wgs84')/1000.0;
  Eclipse_ray_N(rayId).ground_range = ground_range;
end

fprintf('\n')

% plot the rays

for ii = 1:1:num_raybears

  plot3(Eclipse_ray_N(ii).lat,  Eclipse_ray_N(ii).lon, Eclipse_ray_N(ii).height, 'b');
  hold on
end  
legend(UT_str_eclipse, freq_str_eclipse);
hold off
grid on
xlabel('latitude (deg)')
ylabel('longitude (deg)')
zlabel('Height (km)')
cnt=cnt+1;
end



