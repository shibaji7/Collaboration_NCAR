%
% Name:
%   wgs84_xyz2llh_approx.m
%
% Purpose:
%   Convert geocentric x, y, z coordinates to WGS84 ellipsoidal geodetic lat, 
%   long, and height. From the centre of the Earth, x passes through the equator
%   at the prime meridian, y through the equator at longitude of +90 degrees, 
%   and z through the geographic north pole. The conversion routine is an 
%   approximation. Errors are zero for the cases: (1) at all altitudes for 
%   geodetic latitudes of 0 and 90 degrees, and (2) at all latitudes for an 
%   altitude of 0 km. Errors maximise at geodetic latitude of 45 degrees
%   and increases with altitude. For a geodetic latitude of 45 degrees and 
%   altitude of 400km, the error in height is approx. -13cm and in latitude is
%   approx. 1e-8 degrees.
% 
% Calling sequence:
%   [geod_lat, geod_long, height] = wgs84_xyz2llh_approx(x, y, z)
%
% Inputs (scalar):
%   x, y, z    -  geocentric x,y,z (m)
%
% Outputs (scalar):
%   geod_long, geod_lat -  geodetic longitude, latitude (deg)
%   height              -  height above ellipsoid (m)              
% 
% Modification history:
%   24/03/2010 M. A. Cervera
%     Initial version. 
%

function [geod_lat, geod_long, height] = wgs84_xyz2llh_approx(x, y, z)
 
  rtod = 180.0./pi;	       % Radians to degrees.
  a = 6378137.0;               % WGS84 semimajor axis (m)
  f = 1.0 ./ 298.257223563;    % WGS84 flattening factor
  one_minus_f = (1 - f);
  one_minus_f_sq = one_minus_f.^2;

  r = sqrt(x.^2 + y.^2 + z.^2);
  gc_lat = asin(z/r);
  
  xa = one_minus_f*a / sqrt(tan(gc_lat).^2 + one_minus_f_sq);
  mu_a = atan(sqrt(a.^2 - xa.^2) / (one_minus_f*xa));
  if gc_lat < 0, mu_a = -mu_a; end
   
  ra = xa ./ cos(gc_lat);
  l = r - ra;
  
  del_lat = mu_a - gc_lat;
 
  height = l * cos(del_lat);
  
  rho_a = a*one_minus_f_sq / (1 - (2*f - f.^2)*sin(mu_a).^2).^1.5;
  del_mu = atan(l*sin(del_lat)/(rho_a + height));
  
  geod_lat = (mu_a - del_mu) * rtod;
  geod_long = atan2(y, x) * rtod;
  
  return
  
end
