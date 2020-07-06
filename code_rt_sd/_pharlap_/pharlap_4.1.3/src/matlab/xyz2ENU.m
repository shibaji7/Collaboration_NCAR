%
% Name:
%  xyz2ENU.m
%
% Purpose:
%   Convert the location of a point specified in a cartesian coordinate (x,y,z) 
%   frame at a local origin on the Earth to East, North, Up (ENU), where
%   the x axis of the cartesian coordinate frame passes through the equator
%   at the prime meridian, y through the equator at longitude of +90 degrees,
%   and z through the geographic north pole. The (x,y,z) frame origin on the
%   Earth is specified by its geodetic (WGS84) latitude and longitude. 
%
% Calling sequence:
%   [E, N, U] = xyz2ENU(x, y, z, lat, lon)
%
% Inputs: 
%   x, y, z   - cartesian coordinates of the point (m) relative to frame origin,
%               may be arrays of any (but identical) shape
%   lat, lon  - Geodetic (WGS84) latitude and longitude (degrees) of location  
%               of the (x,y,z) frame origin. If x,y,z are arrays then lat, lon
%               must either be be scalar or arrays of identical shape to x,y,z.
%               If x,y,z are scalar then lat, lon may be arrays of any (but
%               identical) shape 
%
% Outputs:
%   E, N, U   - East, North, Up coordinates of point of interest (m) relative
%               to the (x,y,z) frame origin
% 
% Dependencies:
%   None.
%
% Modification history:
%   12/10/2012  V1.0  M. A. Cervera
%     Initial version.
%
%   14/11/2012  V1.1  M. A. Cervera
%     Inputs can now be arrays.
%

function [E, N, U] = xyz2ENU(x, y, z, lat, lon)
 
  % Define constants
  deg2rad = pi/180.0;         % radians to degrees conversion
     
  % determine the sin and cosine of lat and lon required for the series of 
  % rotations to convert local cartesian frame to ENU frame
  lat_r = lat .* deg2rad;
  lon_r = lon .* deg2rad;
  sin_phi = sin(lat_r);
  cos_phi = cos(lat_r);
  sin_theta = sin(lon_r);
  cos_theta = cos(lon_r);

%   % define the rotation matricies 
%   Rot1 = [1  0        0        ; 
%           0  sin_phi  cos_phi  ; 
% 	    0 -cos_phi  sin_phi]; 
% 	
%   Rot2 = [-sin_theta  cos_theta 0  ; 
%           -cos_theta -sin_theta 0  ;
% 	     0          0         1] ;
% 
%   % form the vector of the input coordinates for the rotation matricies to 
%   % operate upon - make sure we record the shape of the original input
%   size_x = size(x);
%   xyz = [x(:) y(:) z(:)]';
%   
%   % Perform rotations to transform local x,y,z cartesian coordinates of the
%   % point to ENU coordinates  
%   ENU = Rot1 * (Rot2 * xyz);
%   
%   % reshape the result to have the same shape as the input 
%   E = squeeze(ENU(1, :));
%   N = squeeze(ENU(2, :));
%   U = squeeze(ENU(3, :));
% 
%   E = reshape(E, size_x);
%   N = reshape(N, size_x);
%   U = reshape(U, size_x);
  
  E = -x.*sin_theta          + y.*cos_theta;
  N = -x.*cos_theta.*sin_phi - y.*sin_theta.*sin_phi + z.*cos_phi;
  U =  x.*cos_theta.*cos_phi + y.*sin_theta.*cos_phi + z.*sin_phi;
  
  
  return  
end
