%
% Name : 
%   ground_fs_loss.m
%
% Purpose :
%   Calculates forward ground scattering losses. Terrain type (land or sea)
%   is taken into account. 
%
% Calling sequence :
%   fs_loss = ground_fs_loss(lat, lon, elev, freq);
%
% Inputs :
%   lat      - latitude of start point (deg)  
%   lon      - longitude of start point (deg)
%   elev     - elevation of the scattered ray (deg)
%   freq     - radio frequency of the ray MHz)
%
% Outputs :
%   fs_loss - the power loss of the scattered radio-waves (dB)
%
% Notes :
%   For propagation purposes the elevation of the scattered ray is equal to
%   the incoming elevation of the ray (i.e. angle of incidence  = angle of
%   reflection) 
%
% Author:
%   V1.0  M.A. Cervera  11/09/2006
%


% This a Matlab help file only. The actual programme is a mex wrapper
% (ground_fs_loss_matlab_wrapper.for) to the Fortran code
% (forward_scatter_loss.for). 
