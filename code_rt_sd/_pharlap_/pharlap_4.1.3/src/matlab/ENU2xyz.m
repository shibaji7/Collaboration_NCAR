%
% Name:
%   ENU2xyz.m
%
% Purpose:
%   Convert the location of a point specified in an East, North, Up (ENU) frame
%   at a local origin on the Earth to cartesian coordinates (x, y, z), where
%   the x axis of the cartesian coordinate frame passes through the equator
%   at the prime meridian, y through the equator at longitude of +90 degrees,
%   and z through the geographic north pole. The ENU frame origin on the Earth 
%   is specified by its geodetic (WGS84) latitude and longitude. 
%
% Calling sequence:
%   [x, y, z] = ENU2xyz(E, N, U, lat, lon)
%
% Inputs: 
%   E, N, U   - East, North, Up coordinates of point of interest (m),  may be
%               arrays of any (but identical) shape 
%   lat, lon  - Geodetic (WGS84) latitude and longitude (degrees) of location  
%               of the (E,N,U) frame origin. If E,N,U are arrays then lat, lon
%               must either be be scalar or arrays of identical shape to E,N,U.
%               If E,N,U are scalar then lat, lon may be arrays of any (but
%               identical) shape 
%
% Outputs:
%   x, y, z  - cartesian coordinates of the point (m) relative to the ENU
%              frame origin
% 
% Dependencies:
%   None.
%
% Modification history:
%   27/11/2009  V1.0  M. A. Cervera
%     Initial version.
%
%   14/11/2012  V1.1  M. A. Cervera
%     Inputs can now be arrays.
%

function [x, y, z] = ENU2xyz(E, N, U, lat, lon)
 
  % Define constants
  deg2rad = pi/180.0;         % radians to degrees conversion
     
  % determine the sin and cosine of lat and lon required for the series of 
  % rotations to convert local ENU frame to cartesian frame
  lat_r = lat .* deg2rad;
  lon_r = lon .* deg2rad;
  sin_phi = sin(lat_r);
  cos_phi = cos(lat_r);
  sin_theta = sin(lon_r);
  cos_theta = cos(lon_r);

  % Perform rotations to tranform local ENU coordinates of the point to
  % x,y,z cartesian coordinates
  x = -E.*sin_theta - N.*sin_phi.*cos_theta + U.*cos_phi.*cos_theta;
  y =  E.*cos_theta - N.*sin_phi.*sin_theta + U.*cos_phi.*sin_theta;
  z =  N.*cos_phi   + U.*sin_phi;

  return  
end
