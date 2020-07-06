%
% Name :
%   irreg_strength.m
%
% Purpose :
%   Simple model of ionospheric irregularity strength in the equatorial and 
%   auroral regions. The returned irregularity strength is the ratio of
%   irregular electron density to the background value at the "phase screen
%   height" of the irregularities.
%
% Calling sequence :
%   strength = irreg_strength(lat, lon, UT, kp)
%
% Inputs (all scalar):
%   lat  - latitude of point (degrees)  
%   lon  - longitude of point (degrees)
%   UT   - 5x1 array containing UTC date and time - year, month, day, 
%          hour, minute
%   kp   - Kp index - used for the auroral irregularity model
%
% Outputs (all scalar):
%   irreg_strength  - ratio of irregular electron density to the background 
%                     value
%
% Note:
%   Currently this routine is a simple model. It is used primarily to 
%   demonstrate the field aligned back scatter capability of
%   raytrace_2d. Treat this model with caution. See M.A. Cervera for further 
%   details. 
%
% Author:
%   V1.0  M.A. Cervera  16/06/2006
%


% This a Matlab help file only. The actual programme is a mex wrapper
% (irreg_strength_matlab_wrapper.f90) to the Fortran code (irreg_strength.f90).
