%
% Name : 
%   dop_spread_eq.m
%
% Purpose :
%   Simple model of Doppler spread imposed on ray traversing the equatorial 
%   region.
%
% Calling sequence :
%   dop_spread = dop_spread_eq(lat, lon, UT, R12);
%
% Inputs :
%   lat      - latitude of start point (deg)  
%   lon      - longitude of start point (deg)
%   UT       - 5x1 array containing UTC date and time - year, month, day, 
%                hour, minute
%   R12      - R12 index
%
% Outputs :
%   dop_spread - square of frequency spread (Hz^2) per unit path length (Km) 
%                  at 1MHz scaled by the electron density (electrons per cm^3) 
%                  squared
%
% Note:
%   Currently this routine is a simple model. It is used primarily to 
%   demonstrate the Doppler spread calculation capability of raytrace_2d.
%   Treat this model with caution. See M.A. Cervera for further details. 
%
% Author:
%   V1.0  M.A. Cervera  16/06/2006
%


% This a Matlab help file only. The actual programme is a mex wrapper
% (dop_spread_matlab_wrapper.for) to the Fortran code (dop_spread_eq.for).
