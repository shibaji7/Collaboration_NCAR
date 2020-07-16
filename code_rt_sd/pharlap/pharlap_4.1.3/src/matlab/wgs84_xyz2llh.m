%
% Name:
%   wgs84_xyz2llh.m
%
% Purpose:
%   Convert geocentric x, y, z coordinates to WGS84 ellipsoidal geodetic lat, 
%   long, and height.
% 
% Calling sequence:
%   [geod_lat, geod_long, height] = wgs84_xyz2llh(x, y, z)
%
% Inputs (scalar):
%   x, y, z    -  geocentric x,y,z (m)
%
% Outputs (scalar):
%   geod_long, geod_lat -  geodetic longitude, latitude (deg)
%   height              -  height above ellipsoid (m)              
% 
% Modification history:
%   31/07/2007 M. A. Cervera
%     Initial version. Converted from IDL code written by  R. Sterner,   
%     Johns Hopkins University/Applied Physics Laboratory, 2002 May 06. 
%

function [geod_lat, geod_long, height] = wgs84_xyz2llh(x, y, z)
 
  rtod = 180.0./pi;	       % Radians to degrees.
  a = 6378137.0;               % WGS84 semimajor axis (m)
  f = 1.0 ./ 298.257223563;    % WGS84 flattening factor
  b = a.*(1-f);                % WGS84 semiminor axis (m)

  a2 = a.^2;
  a4 = a.^4;
  b2 = b.^2;
  b4 = b.^4;
  a2b2 = a2.*b2;

  % Convert problem from 3-D to 2-D
  x0 = sqrt(x.^2 + y.^2);
  y0 = z;
  x02 = x0.^2;
  y02 = y0.^2;

  % Coefficients of the polynomial
  c = zeros(1, 5);
  c(5) = a2b2.*(b2.*x02 + a2.*y02 - a2b2);
  c(4) = 2.*a2b2.*(x02 + y02 - a2 - b2);
  c(3) = a2.*(x02 - 4.*b2) + b2.*y02 - a4 - b4;
  c(2) = -2.*(a2 + b2);
  c(1) = -1.0d0;

  % Find roots of the 4th degree polynomial.   
  rr = roots(c);

  % Nearest = nadir point.
  t = rr(4);	

  % Ellipse X, Y at nadir point.
  xe = a2.*x0 ./ (a2 + t);	   
  ye = b2.*y0 ./ (b2 + t);

  % Calculate Geodetic height (m)
  height = sign(t) .* sqrt((x0 - xe).^2 + (y0 - ye).^2);

  % Calculate geocentric latitude (radians).
  lat0r = atan2(ye, xe);

  % Calculate Geodetic lat and long
  geod_lat_r = atan(tan(lat0r) .* a2 ./ b2);
  geod_long_r = atan2(y, x);		
  geod_lat = geod_lat_r .* rtod;	        % Lat in degrees
  geod_long = geod_long_r .* rtod;		% Long in degrees

  return
  
end
