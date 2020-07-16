%
% Name : 
%   ground_bs_loss.m
%
% Purpose :
%   Calculate power loss (dB) of the radio-waves back-scattered from the ground.
%   This is a simple model of back-scattered loss which is based on two land
%   types viz. land or sea. For the case of sea, the sea is considered to be
%   fully-developed and the back scatter loss is 26dB. For the case of land
%   the back scatter loss is 29 dB.
%
% Calling sequence :
%   bs_loss = ground_bs_loss(lat, lon);
%
% Inputs :
%   lat      - latitude of start point (deg)  
%   lon      - longitude of start point (deg)
%
% Outputs :
%   bs_loss - the power loss of the radio-waves back-scattered from the 
%             ground (dB)
%
% Author:
%   V1.0  M.A. Cervera  11/09/2006
%


% This a Matlab help file only. The actual programme is a mex wrapper
% (ground_bs_loss_matlab_wrapper.for) to the Fortran code (land_type.for).
