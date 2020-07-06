Refer script pharlap.m for PHaRLAP toolbox documentation and details of all the subroutines that can be used.
The following scripts for used the purpose of my project that are to be run from pharlap_4.1.3 folder:
For 2d raytracing:
gen_iono_grid_2d_sneh.m
raytrace_2d.m
Modified scripts:
Eclipse_raytracing.m, NoEclipse_raytracing.m
New additions: 
runThis.m: To specify inputs to functions Eclipse_raytracing.m, NoEclipse_raytracing.m

For 3d raytracing:
gen_iono_grid_3d modified to gen_iono_grid_3d_new
raytrace_3d.m
raz2latlon.m: converts range and azimuth from origin to spherical earth latitude and longitude for various geoids
New additions: 
Spline_3d.m: Generates a 3d attenuation function in latitude and longitude and computation to scale latitude and and longitude axis to the corresponding ionospheric subgrid (eclipsed region) 
Set_ionogrid_3d.m: To obtain 3d ionopheric grid parameters for uneclipsed and eclipsed cases
NoEclipse_3d_azim_new.m: 3d raytrace plots for uneclipsed, normal ionosphere
Eclipse_3d_azim_new.m: 3d raytrace plots for eclipsed ionosphere
plot_azim_3d.m: Uneclipsed and eclipsed ray traces for varying azimuths and fixed elevation on a single plot for comparison 
plot_elev_3d.m: Uneclipsed and eclipsed ray traces for varying elevations and fixed azimuth on a single plot for comparison 
runThis_3d.m: 3d raytracing in latitude, longitude and height for eclipsed and uneclipsed ionosphere
