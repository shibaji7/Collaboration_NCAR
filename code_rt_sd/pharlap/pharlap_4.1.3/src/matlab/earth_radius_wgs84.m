%
% Name:
%   earth_radius_wgs84.m
%
% Purpose:
%   Returns the "radius" of the Earth (i.e. distance from centre of the Earth
%   to the surface) at a given input WGS84 geodetic latitude.
%	
% Calling sequence:
%     re_wgs84 = earth_radius_wgs84(geod_lat)
%
% Inputs:	
%   geod_lat - WGS84 geodetic latitude (degrees)
%
% Outputs:
%   re_wgs84 - "radius" of the Earth at the input geodetic latitude (m)
%
% Dependencies:
%   none
%
% Modification History:
% 08/11/2007 M. A. Cervera
%   Initial version. 
%

function re_wgs84 = earth_radius_wgs84(geod_lat)

  sem_maj_ax = 6378137.0;         % semi-major axis of WGS84 (m)
  sem_min_ax = 6356752.314245;    % semi-minor axis of WGS84 (m)
  dtor = pi / 180.0;

  a2 = sem_maj_ax.^2;
  b2 = sem_min_ax.^2;
  a2cos2l  = a2 * cos(geod_lat * dtor).^2;
  b2sin2l  = b2 * sin(geod_lat * dtor).^2;

  re_wgs84 = sqrt((a2 * a2cos2l + b2 * b2sin2l) / (a2cos2l + b2sin2l));

end
