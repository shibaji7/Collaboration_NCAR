% 
% Name :
%   coning.m
% 
% Purpose :
%   Calculates correction to ray azimuth due to cone effect of linear
%   arrays. If the radar antenna beam direction is given by "radar_bearing" and
%   antenna array bore direction is "radar_bore", then the "ray_bearing will" be
%   given by 
%
%       ray_bearing = radar_bearing + coning_correction 
%
%   where coning_correction is the correction due to the cone effect and given by
%
%       coning_correction = coning(elev, radar_bearing - radar_bore)
%    
%   and "elev" is the ray elevation.
%
% Calling sequence :
%    coning_correction = coning(elev, off_bore);
%
% Inputs :
%   elev - the ray elevation (degrees, scaler)
%   off  - the azimuth of the radar bearing from radar bore (degrees, scaler)
%
% Outputs :
%   coning_correction - the correction to the ray azimuth (degrees)
% 
% Modification History :
% 26/10/2005  V1.0  M. A. Cervera  Author.
%
% 07/04/2009  V1.1  M. A. Cervera
%   Added error checking to make sure that the input elevation does not
%   exceed 90 - off_bore angle.
%

function coning_correction = coning(elev, off_bore)

deg2rad = pi ./ 180;   % degrees to radians conversion
rad2deg = 180 ./ pi;   % radians to degrees conversion

if (elev > 90-off_bore) error('input elevation is > 90 - off_bore angle'), end
  
sin_elev = sin(deg2rad .* elev);
sin_off = sin(deg2rad .* off_bore);
temp = sqrt(abs(1.e0 - sin_off.^2 - sin_elev.^2));
coned = 90 - atan2(temp, sin_off) .* rad2deg;
coning_correction = coned - off_bore;

return

