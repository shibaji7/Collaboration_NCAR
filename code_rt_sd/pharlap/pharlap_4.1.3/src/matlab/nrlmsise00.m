%
% Name :
%   nrlmsise00.m
%
% Purpose :
%   Matlab wrapper to the NRLMSISe-00 model atmosphere routines distributed
%   with the IRI-2012 Fortan based empirical model ionosphere. 
%   Returns modeled atmospheric temperature and densities.
%
% Calling sequence :
%    [d, t] = nrlmsise00(lat, lon, height, UT)
%    [d, t] = nrlmsise00(lat, lon, height, UT, R12)
%
% Inputs :
%   lat      - geographic latitude of point (degrees)  
%   lon      - geographic longitude of point (degrees)
%   height   - height (km)
%   UT       - 5x1 array containing UTC date and time - year, month, day, 
%                hour, minute
% Optional Inputs:
%   R12      - yearly-smoothed monthly median sunspot number. If omitted or
%              set to -1 the R12 value will be read from file (ig_rz.dat_
%              and my be historical or projected conditions (dependent on
%              the epoch).
%
% Outputs :
%    d - 9xN array of densities, all in cm-3 except for d(6)
%      d(1, :) = He number density
%      d(2, :) = O number density
%      d(3, :) = N2 number density
%      d(4, :) = O2 number density
%      d(5, :) = Ar number density
%      d(6, :) = Total mass density (g/cm3)
%      d(7, :) = H number density
%      d(8, :) = N number density
%      d(9, :) = Anomalous oxygen number density
%    t - 2xN array of temperatures in K
%      t(1, :) = Exospheric temperature
%      t(2, :) = Temperature at specified height
%
% Author:
%   V1.0  L.H. Pederick  30/08/2012
%


% This a Matlab help file only. The actual programme is a mex wrapper
% (nrlmsise00_matlab_wrapper.c) to the Fortran code (cira.for).
