%
% Name :
%   solar_za.m
%
% Purpose :
%   Calculate the solar zenith angle (apprximate formula).
%
% Calling sequence :
%   solar_zen_ang = solar_za(lat, lon, UT)
%
% Inputs :
%   lat - latitude (degrees)
%   lon - longitude (degrees)
%   UT  - 5xN array containing UTC date and time - year,month,day,hour,minute
%
% Outputs :
%   solar_zen_ang - solar zenith angle, degrees
%
% Author:
%   V1.0  M.A. Cervera  07/08/2009
%
%   V1.1  L.H. Pederick 18/09/2014
%      Modified to handle vectorized input for UT

function solar_zen_ang = solar_za(lat, lon, UT)
  
  if size(UT,1) == 1, UT = UT'; end

  % convert lat, lon to radians
  dtor = pi ./ 180;
  lat_r = lat .* dtor;
  lon_r = lon .* dtor;
  
  % calculate the hour angle, hour_ang
  hour = UT(4,:) + UT(5,:)./60;                   % UTC decimal hour
  hour_LST = mod(hour + lon./15, 24);         % Local Solar Time decimal hour 
  hour_ang = 15 .* (hour_LST - 12) .* dtor;   % Hour angle
  
  % calculate the day number in the year, doy
  doy = julday(UT(3,:), UT(2,:), UT(1,:)) - julday(zeros(1, size(UT,2)), ones(1, size(UT,2)), UT(1,:));
  
  % calculate the solar declination, solar_dec (approximate formula)
  obliq_ecliptic = 23.45;
  solar_dec = obliq_ecliptic .* sin(2.*pi .* (doy + 284) ./ 365) .* dtor; 
  
  % calculate the solar zenith angle in degrees
  csza = sin(lat_r) .* sin(solar_dec) + ...
         cos(lat_r) .* cos(solar_dec) .* cos(hour_ang);
  solar_zen_ang = acos(csza) ./ dtor;

  return
end

