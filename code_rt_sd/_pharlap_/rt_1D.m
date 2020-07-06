%% Description
%
% Author : Chakraborty, Shibaji
% Name : 1D Ray trace model for SWF
%
% Purpose :
%   Estimate the SuperDARN radar rays (for only one beam) due to
%   change in ionospheric conditions during solar flare events using the
%   the spherical Earth 1D NRT raytrace_2d_sp for a fan of rays.
%   Ray trajectories are saved under 'local' directory
%
% Calling sequence :
%   rt_1d
%
% Inputs :
%   Radar code
%   Time index
%   Beam num
%   Directory
%   
% Outputs :
%   .mat files with ray trajectory


%% Close and clear all panes and variables and command line window and load PharLap
close all
startup;

%% initialize
R12 = 100;
speed_of_light = 2.99792458e8;
doppler_flag = 1;
irregs_flag = 0;
kp = 0;

%% load data
load(['../' dic 'bearing.mat'])
load(['../' dic 'tmp.mat'])

%% initialize IRI grid
plat
clear iri_options
iri_options.Ne_B0B1_model = 'Bil-2000';
tic
%[iono_pf_grid, iono_pf_grid_5, collision_freq, irreg] = ...
%			gen_iono_grid_2d(olat, olon, R12, UT, rb, ...
%			max_range, num_range, range_inc, start_height, ...
%			height_inc, num_heights, kp, doppler_flag, 'iri2016', ...
%			iri_options);
toc
