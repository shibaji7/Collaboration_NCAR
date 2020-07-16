%
% Name:
%   wgs842gc_lat.m
%
% Purpose:
%   Convert WGS84 ellipsoidal geodetic latitude to geocentric latitude at a 
%   given geodetic height above the ground. 
%      NB 1. Geocentric longitude and geodetic longitude have the same value. 
%         2. Geodetic latitude is defined wrt the normal to the ellipsoid and
%            is constant with height above the ellipsoid.
%         3. Geocentric latitude is defined wrt the centre of the Earth. Thus it 
%            NOT cons.tant with height above the ellipsoid (only true for a 
%            spheroid)
% 
% Calling sequence:
%   gc_lat = wgs842gc_lat(geod_lat, height)
%
% Inputs:
%   geod_lat  -  geodetic longitude (degrees)
%   height    -  geodetic height above the ground (m)
%
% Outputs:
%   gc_lat  -  geocentric latitude (degrees)
% 
% Dependencies:
%   none
%
% Modification history:
%   08/11/2007 M. A. Cervera
%     Initial version.
%

function gc_lat = wgs842gc_lat(geod_lat, height)

  rtod = 180.0d0/pi;
  dtor = pi/180.0d0;
  a = 6378137.0d0;              % WGS84 semimajor axis (m)
  f = 1.0 / 298.257223563d0;    % WGS84 flattening factor
  b = a*(1-f);	                % WGS84 semiminor axis (m)
  e2 = (a.^2 - b.^2) / a.^2;    % eccentricity of ellipsoid squared
     
  % Convert input WGS84 geodetic lat, long to radians
  geod_lat_r = geod_lat * dtor;

  % Calculate the geocentric latitude
  chi = sqrt(1 - e2 * sin(geod_lat_r).^2);
  c1 = a + chi * height;
  tan_gc_lat = tan(geod_lat_r) * (c1 - a * e2) / c1;
  gc_lat = atan(tan_gc_lat)* rtod;

end
