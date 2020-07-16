% 2d raytracing for west, north and east rays starting from Fort Hays location
clc;
clear all
clearvars;
close all;
% Fort Hays
origin_lat_p = '38.86';           % latitude of the start point of ray
origin_long_p = '-99.39';         % longitude of the start point of ray
% North
ray_bear_p = '0';
UT_Array_n = [2016 8 21 17 58]; % Universaltime of the centre of the eclipse along the direction of the ray 
dir_path = 'raytrace-north';
offset_n = 240;                 % Distance of the centre of the eclipse from origin (Fort Hays) along the direction of the ray 
NoEclipse_Raytracing(origin_lat_p,origin_long_p,ray_bear_p,UT_Array_n, dir_path);
Eclipse_Raytracing(origin_lat_p,origin_long_p,ray_bear_p, UT_Array_n, dir_path,offset_n);
close all;
% West
ray_bear_p = '-69'; 
UT_Array_n = [2016 8 21 17 37];
dir_path = 'raytrace-west';
offset_n = 958;
NoEclipse_Raytracing(origin_lat_p,origin_long_p,ray_bear_p, UT_Array_n, dir_path);
Eclipse_Raytracing(origin_lat_p,origin_long_p,ray_bear_p,UT_Array_n, dir_path,offset_n);
close all;
% East
ray_bear_p = '89';
UT_Array_n = [2016 8 21 18 14];
dir_path = 'raytrace-east';
offset_n = 534;
NoEclipse_Raytracing(origin_lat_p,origin_long_p,ray_bear_p,UT_Array_n, dir_path);
Eclipse_Raytracing(origin_lat_p,origin_long_p,ray_bear_p, UT_Array_n,dir_path,offset_n);
%Movie(10);