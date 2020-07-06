%
% Name : 
%   gm_freq_offset.m
%
% Purpose :
%   Calculates the approximate geomagnetic O-X mode frequency split (MHz) for
%   a specified propagation path. This can then be applied to a ray
%   calculated using ART, 2D NRT or "no geomagnetic field" 3D NRT to give an
%   approximate O and X mode solution.
%
% Calling sequence :
%   [del_freq_O, del_freq_X] = ...
%      gm_freq_offset(lat_cp, lon_cp, ht_cp, ray_bearing, freq, plas_freq, UT)
%
% Inputs :
%   lat_cp      - latitude of control point of ray (i.e. the ray apogee) (deg)  
%   lon_cp      - longitude of control point of ray (deg)
%   ht_cp       - height of control point of ray (deg)
%   ray_bearing - bearing of ray from tx. (deg)
%   freq        - radio frequency of the ray (MHz)
%   UT          - universal time, 5x1 array (year, month, day, hour, minute)
%   plas_freq   - plasma freq of the ionosphere at the ray control point (MHz)
%
% Outputs :
%   del_freq_O  -  shift to be applied to ray frequency for O mode case (MHz)
%   del_freq_X  -  shift to be applied to ray frequency for X mode case (MHz)
%
% Notes :
%   1. The "no geomagnetic field" ray must also have the group path offset
%      calculated. If group_path, phase_path and freq are the "no field"
%      values, the the group_path offsets are:
%
%      del_group_O = (group_path - phase_path) .* del_freq_O / freq;
%      del_group_X = (group_path - phase_path) .* del_freq_X / freq;
%
%   2. Application of the calculated group path and frequency shifts is as
%      follows :
%
%      group_O = group + del_group_O;
%      group_X = group + del_group_X;
%
%      freq_O = freq + del_freq_O;
%      freq_X = freq + del_freq_X;
%
%   3. see the example code ois_synth.m for an example on how gm_freq_offset.m
%      is used
%
%
% References :
%   1. Bennett et al. (1994), "Analytic calculation of the ordinary (O) and 
%      extraordinary (X) mode nose frequencies on oblique ionograms", J. Atmos.
%      Terr. Phys., Vol 56, No. 5, pp 631-636.
%
%   2. Dyson and Bennett (1980), "A Universal Chart for Use in Relating
%      Ionospheric Absorption to Phase Path and Group Path", IEEE
%      Transactions on Antennas and Propagation, Vol AP-28, No. 3, 380-384.
%
%   3. Bennett et al. (1991), "Analytic Ray Tracing for the Study of HF
%      Magneto-ionic Radio Propagation in the ionosphere", Applied
%      Computational Electromagnetics Society Journal, Vol 6(1), 192-210
%      (see Appendix)
%
%   4. Chen et al. (1990), "Automatic fitting of quasi-parabolic segments to
%      ionospheric profiles", J. Atmos. Terr. Phys, Vol 52, No. 4, 277-288 
%      (see Appendix B)
%
%
% Modification history:
%   V1.0  M.A. Cervera  08/05/2015
%       Based on Matlab code and algorithm provided by Andrew Heitman
%

function [del_freq_O, del_freq_X] = ...
      gm_freq_offset(lat_cp, lon_cp, ht_cp, ray_bearing, freq, plas_freq, UT)

			      

% check inputs
if (~isscalar(lat_cp) | ~isscalar(lon_cp)  | ~isscalar(ht_cp)  | ...
    ~isscalar(ray_bearing) | ~isscalar(freq) | ~isscalar(plas_freq))
  error('Inputs 1 - 6 must be scalars')
  return
end

% Return frequency offsets of zero when input frequency is zero.
if (freq == 0)
  del_freq_O = 0;
  del_freq_X = 0;
  
  return
end

% Some useful constants for easy conversion between degrees and radians.
dtor = pi/180;
rtod = 180/pi;

% gyro_factor - converts magnetic field strength (Tesla) to electron
% gyro-frequency 
gyro_factor = 2.7992492071667d4; 

%  geomagnetic field parameters at the control point.
mag_field = igrf2016(lat_cp, lon_cp, UT, ht_cp);
B_mag = mag_field(4);
B_dip = mag_field(8);
B_dec = mag_field(10);

gyro_freq = B_mag .* gyro_factor;

% angle between B-field and ray direction at control point
theta = acos(cos(B_dip .* dtor) .* cos((B_dec - ray_bearing) .* dtor));


% Compute o-mode and x-mode offsets from the effective frequency, as defined
% by Bennett, Chen & Dyson (1991). Strictly speaking, 'freq' is assumed to
% represent this unperturbed value.

% X = (fN/f)^2 and Y = fH/f are the normalised frequency variables
% W = the ratio of X-1 and Y (negative for HF ionospheric propagation)
Y = gyro_freq ./ freq;
X = (plas_freq ./ freq).^2;
W_o = (X - 1) ./ Y;
W_x = (X - sqrt(X).*Y - 1) ./ Y;

% Calculate o-mode and x-mode h factors at the ray's control point. 
% Theta is the angle between the wave normal (i.e. direction of phase
% propagation) and the geomagnetic field. As we can't calculate this value 
% (full 3D magneto-ionic raytrace is required) we will approximate this value 
% with the angle between the ray direction and the geomagnetic field. This
% approximation is OK for oblique propagation but errors will grow as the ray
% elevation -> 90 deg. 
h_o = h_param(W_o, theta, +1);
h_x = h_param(W_x, theta, -1);

% Calculate the effective frequency offsets.
g_o = Y .* h_o/2;
g_x = Y .* h_x/2;
del_freq_O = freq .* (g_o ./ (g_o - 1));  
del_freq_X = freq .* (g_x ./ (g_x - 1)); 


return
end


%-------------------------------------------------------------------------------

function h = h_param(W, theta, OX_mode)

% Equation 7 from Dyson & Bennett (1980).

% The input argument OX_mode specifies the O (+1) or X mode (-1) version of
% the equation. 

% stay away from theta = 90 deg as denominator -> 0 here for X mode
idx = find(theta > pi/2 - 1e-5 & theta < pi/2 + 1e-5);
theta(idx) = theta(idx) - 1e-5;

S2 = sin(theta).^2;
C2 = cos(theta).^2;

fact = sqrt(1 + 4*C2 .* W.^2 ./ S2.^2);
numerator = 2*C2 .* W .* (1 - OX_mode*fact);
denominator = 1 + C2 + OX_mode*S2.*fact + 4*C2 .* W.^2 ./ S2;
h = numerator ./ denominator;

% at the limits theta -> 0 and theta -> 180 :
%    h -> +1 for O mode if W != 0
%          0 for O mode if W == 0
%    h -> -1 for X mode for all W
idx = find(S2 < 1e-10);
if OX_mode == -1
  h(idx) = OX_mode;   
else
  if (W == 0)
    h(idx) = 0;
  else
    h(idx) = 1;
  end
end

return
end