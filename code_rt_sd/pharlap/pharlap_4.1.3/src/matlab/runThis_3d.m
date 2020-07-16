% Script to obtain 3d raytracing in latitude, longitude and height for eclipsed and uneclipsed ionosphere
% First generates the 3d attenuation factors to be used for the eclipsed case
% Output: Azimuth plots

clc;
clear all
clearvars;
close all;
run('spline_3d.m');
run('set_ionogrid_3d.m');