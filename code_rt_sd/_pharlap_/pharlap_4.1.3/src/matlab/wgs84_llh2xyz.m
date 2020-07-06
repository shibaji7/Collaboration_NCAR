%
% Name:
%   wgs84_llh2xyz.m
%
% Purpose:
%   Convert WGS84 ellipsoidal geodetic lat, long, and height to geocentric 
%   x, y, z coordinates.
% 
% Calling sequence:
%   [x, y, z] = wgs84_llh2xyz(geod_lat, geod_long, height)
%
% Inputs:
%   geod_lat  -  array of geodetic latitude (deg)
%   geod_long -  array of geodetic longitude
%   height    -  array of heights above ellipsoid (m)              
% 
% Outputs: 
%   x, y, z    -  arrays of geocentric x, y, z (m)
%
% Dependencies:
%   none
%
% Modification history:
%   01/08/2007 M. A. Cervera
%     Initial version.  
%

function [x, y, z] = wgs84_llh2xyz(geod_lat, geod_long, height)

  dtor = pi./180.0;	        % Degrees to radians.
  a = 6378137.0;                % WGS84 semimajor axis (m)
  f = 1.0 ./ 298.257223563;     % WGS84 flattening factor
  b = a.*(1-f);                 % WGS84 semiminor axis (m)
  e2 = (a.^2 - b.^2) ./ a.^2;   % eccentricity of ellipsoid squared

  % Convert geodetic lat, long to radians
  geod_long_r = geod_long .* dtor;		
  geod_lat_r = geod_lat .* dtor;
  
  % do the calculation
  sin_gd_lat = sin(geod_lat_r);
  cos_gd_lat = cos(geod_lat_r);
  sin_gd_lon = sin(geod_long_r);
  cos_gd_lon = cos(geod_long_r);

  chi = sqrt(1 - e2 .* sin_gd_lat.^2);

  x = (a./chi + height) .* cos_gd_lat .* cos_gd_lon;
  y = (a./chi + height) .* cos_gd_lat .* sin_gd_lon;
  z = ((a .* (1 - e2)) ./ chi + height) .* sin_gd_lat;

  return
  
end
