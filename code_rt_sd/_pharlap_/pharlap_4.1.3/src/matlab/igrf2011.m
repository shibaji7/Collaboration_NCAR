%
% Name :
%   igrf2011.m
%
% Purpose :
%   Matlab wrapper to the IGRF geomagnetic field routines (igrf.for)
%   distributed with the IRI-2011 fortan based empirical model ionosphere. 
%   Returns modeled geomagnetic field.
%
% Calling sequence :
%    mag_field = igrf2011(lat, lon, UT, height)
%
% Inputs :
%   lat      - geographic latitude of point (degrees)  
%   lon      - geographic longitude of point (degrees)
%   UT       - 5x1 array containing UTC date and time - year, month, day, 
%                hour, minute
%   height  -  height (km)
%
% Outputs :
%   mag_field - 1 X 10 array of IGRF magnetic field parameters
%      mag_field(1) = North component of magnetic field (Tesla) 
%      mag_field(2) = East  component of magnetic field (Tesla) 
%      mag_field(3) = Downwards component of magnetic field (Tesla)
%      mag_field(4) = magnetic field strength (Tesla)
%      mag_field(5) = dipole moment (Tesla)
%      mag_field(6) = L value
%      mag_field(7) = L flag (1 = L is correct; 2 = L is not correct;
%                             3 = approximation is used)
%      mag_field(8) = magnetic dip (degrees)
%      mag_field(9) = dip latitude (or magnetic latitude i.e. atan(tan(dip)/2)
%                     (degrees)
%      mag_field(10) = magnetic declination (degrees)
%
% Author:
%   V1.0  M.A. Cervera  05/09/2012
%


% This a Matlab help file only. The actual programme is a mex wrapper
% (igrf2011_matlab_wrapper.for) to the Fortran code (igrf.for).
