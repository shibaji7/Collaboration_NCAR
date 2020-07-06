% 
% Name : 
%   raytrace_3d_sp
%
% Purpose :
%   3D magneto-ionic numerical raytrace for a multihop ray The computational
%   routine is multi-threaded and accepts a user specified number of independant
%   rays which it farms out to the available computational cores.
%
%   The Hamiltonian for the ray is that set down by Haselgrove and Hasegrove
%   (1960), and Hasegrove (1963). Geomagnetic field effects are considered
%   for O and X  polarized radio waves. No-field case is also considered if
%   required.  Spherical Earth coordinate system is assumed with user input
%   radius of the Earth. This version takes approximately half the time to
%   complete that the WGS84 compliant raytrace_3d does. An appropriate choice
%   for   the radius of the Earth will minimise the error wrt the WGS84
%   coordinate  system (see also earth_radius_wgs84.m)
%
% Calling sequence :
%   The first time raytrace_3d_sp is called, ionospheric and geomagnetic grids
%   must be specified. For subsequent calls if these grids are not passed in 
%   (see calling sequence examples 4 5 and 6) then the grids from the
%   previous call (which are in memory) will be used. This may speed up the 
%   subsequent calls to raytrace for large grids. 
%
%   1.  [ray_data, ray_path_data, ray_state_vec] = ...
%           raytrace_3d_sp(origin_lat, origin_long, origin_ht, elev,ray_bearing,
%                       freq, OX_mode, nhops, tol, rad_earth, iono_en_grid, ...
%    	                iono_en_grid_5, collision_freq, iono_grid_parms, ...
%                       Bx, By, Bz, geomag_grid_parms)
%
%   2.  [ray_data, ray_path_data, ray_state_vec] = ...
%           raytrace_3d_sp(origin_lat, origin_long, origin_ht, elev,ray_bearing,
%                       freq, OX_mode, nhops, tol, rad_earth, iono_en_grid, ...
%    	                iono_en_grid_5, collision_freq, iono_grid_parms, ...
%                       Bx, By, Bz, geomag_grid_parms, ray_state_vec_in)
%
%   3.  [ray_data, ray_path_data, ray_state_vec] = ...
%           raytrace_3d_sp(origin_lat, origin_long, origin_ht, elev, ray_bearing,
%                       freq, OX_mode, nhops, tol, rad_earth)
%
%   4.  [ray_data, ray_path_data, ray_state_vec] = ...
%           raytrace_3d_sp(origin_lat, origin_long, origin_ht, elev,ray_bearing,
%                       freq, OX_mode, nhops, tol, rad_earth, ray_state_vec_in)
%
% Calling sequence (2) : THESE ARE FOR BACKWARDS COMPATIBILITY WITH PHARLAP 
%                        VERSION 3.7.1 AND EARLIER
%   For these calling sequences only a single ray is being traced (i.e. 
%   length(elevs) = length(freqs) = length(ray_bearings) = 1)
%
%   5.  [ray_data,ray_path_data,nhops_attempted,ray_label,ray_state_vec] = ...
%           raytrace_3d_sp(origin_lat, origin_long, origin_ht, elev,ray_bearing,
%                       freq, OX_mode, nhops, tol, rad_earth, iono_en_grid, ...
%    	                iono_en_grid_5, collision_freq, iono_grid_parms, ...
%                       Bx, By, Bz, geomag_grid_parms)
%
%   6.  [ray_data,ray_path_data,nhops_attempted,ray_label,ray_state_vec] = ...
%           raytrace_3d_sp(origin_lat, origin_long, origin_ht, elev,ray_bearing,
%                       freq, OX_mode, nhops, tol, rad_earth, iono_en_grid, ...
%    	                iono_en_grid_5, collision_freq, iono_grid_parms, ...
%                       Bx, By, Bz, geomag_grid_parms, ray_state_vec_in)
%
%   7.  [ray_data,ray_path_data,nhops_attempted,ray_label,ray_state_vec] = ...
%           raytrace_3d_sp(origin_lat, origin_long, origin_ht, elev,ray_bearing,
%                       freq, OX_mode, nhops, tol, rad_earth)
%
%   8.  [ray_data,ray_path_data,nhops_attempted,ray_label,ray_state_vec] = ...
%           raytrace_3d_sp(origin_lat, origin_long, origin_ht, elev,ray_bearing,
%                       freq, OX_mode, nhops, tol, rad_earth, ray_state_vec_in)
%
% Inputs :
%   origin_lat      - geocentric latitude (-90 to 90 degrees) of start point
%                     of rays (All rays have the same origin)
%   origin_long     - geocentric longitude (-180 to 180 degrees) of start
%                     point of rays
%   origin_height   - height above sea-level of start point of ray which
%                     must be below start of ionosphere (km). (If the start of 
%                     the ray tracing is inside the ionosphere then the initial
%                     state vector of the ray must be specified.)
%   elevs           - 1 X M (where M is the number of rays) array of initial
%                     elevation of rays (deg) 
%   ray_bearings    - 1 X M array of initial bearing of the rays (deg) 
%   freqs           - 1 X M array of wave frequency of the rays (MHz)
%   OX_mode         - polarization mode of ray: 1 = O, -1 = X, 0 = no field
%   nhops           - number of hops to complete   
%   tol             - a 1 or 3 element vector controlling ODE solver precision.
%      tol(1)  =  ODE solver tolerence, valid values 1e-12 to 1e-2
%      tol(2)  =  ODE solver minimum step size to consider (0.001 to 1 km)
%      tol(3)  =  ODE solver maximum step size to consider (1 to 100 km)
%      
%      If tol is a scalar then min and max step sizes are set to 0.01 and
%      10km. Suggested value for tol in this case is 1e-7. This provides
%      backward compatibility with PHaRLAP version 3.2.1 and earlier.
%
%      Example values: 1. Highest precision, slowest  - tol(1) = 1e-8
%                                                       tol(2) = 0.01 km
%                                                       tol(2) = 10 km
%                      2. Lower precision, faster     - tol(1) = 1e-7
%                                                       tol(2) = 0.025 km
%                                                       tol(2) = 25 km 
%                      3. Lowest precision, fastest   - tol(1) = 1e-6
%                                                       tol(2) = 0.1 km
%                                                       tol(2) = 100 km 
%    
%      The example values may be selected by setting tol to be a scalar with
%      a value of 1, 2, or 3.
%
%   rad_earth       - radius of Earth to use for coord conv. routines (m)
%   iono_en_grid    - 3D grids of electron density (num/cm^3), elec. density
%   iono_en_grid_5    5 minutes later, and collision frequency (MHz) as a 
%   collision_freq    function of  latitude, longitude and height. The 
%                     maximum grid sizes are 701 X 701 X 301 (lat, lon, height).
%                     (501 x 501 x 301 on 32-bit Windows platforms). The 
%                     ionospheric data and their gradients must be smooth and 
%                     continuous. See note 2 below on memory.
%
%   iono_grid_parms - 9x1 vector containing the parameters which define the
%                     ionospheric grid :
%           (1) geocentric latitude (degrees) of start of grid - must be in the
%               range -90 to 90 degrees 
%           (2) latitude step (degrees)
%           (3) number of latitudes
%           (4) geocentric longitude (degrees) of start of grid - must be in
%               the range -180 to 180 degrees
%           (5) lonfitude step (degrees)
%           (6) number of longitudes
%           (7) geocentric height (km) start
%           (8) height step (km)
%           (9) number of heights
%
%   Bx, By, Bz        - 3D grids of x, y and z components of the geomagnetic 
%                       field (Tesla) as a function of latitude,
%                       longitude and height. The maximum grid size is 
%                       101 X 101 X 201 (lat, lon, height). The geomagetic field
%                       data and their gradients must be smooth and continuous.
% 
%   geomag_grid_parms - 9x1 vector containing the parameters which define the
%                       ionospheric grid :
%           (1) geocentric latitude (degrees) of start of grid - must be in the
%               range -90 to 90 degrees 
%           (2) latitude step (degrees)
%           (3) number of latitudes
%           (4) geocentric longitude (degrees) of start of grid - must be in
%               the range -180 to 180 degrees
%           (5) longitude step (degrees)
%           (6) number of longitudes
%           (7) geocentric height (km) start
%           (8) height step (km)
%           (9) number of heights
%
% Optional Inputs :
%   ray_state_vec_in(:) - 1 X M structure containing the state vector of each
%            ray (total of M rays) at its starting point. Each field (total of
%            11 fields) is an element of the state vector. If the first field
%            for a given ray is -1 then default starting values are used for
%            that ray. NB: non-default values are for EXPERTS ONLY - you really
%            want to know what you are doing! If used then the input starting
%            point, elevation and bearing of ray (first 5 inputs) will be
%            ignored. The fields of this struture are:
%     .pos_x          - cartesian x position of ray start point (m)
%     .pos_y          - cartesian y position of ray start point (m)
%     .pos_z          - cartesian z position of ray start point (m)
%     .dir_x          - cartesian x, y, and z components of vector with a
%                       magnitude equal to the refractive index pointing in the
%                       direction of the wave normal (which, generally, is not 
%                       in the direction  of the ray propagation) 
%     .group_path     - group path (m)
%     .geome_path     - geometrical distance travelled by rays (m)
%     .phase_path     - phase path (m)
%     .indep_var      - independant variable 
%     .ODE_step_size  - ODE solver independant variable step size
%                                                                            
% Outputs :
%   ray_data(:) -  1 X M structure containing the information for each of the
%                  rays (total of M rays). Each field is a 1 X N array
%                  containing information for each hop. The fields of this
%                  structure are: 
%     .lat                   - geocentric latitude of end point of ray (deg) 
%                              for each hop 
%     .lon                   - geocentric longitude of end point of ray (deg)
%                              for each hop 
%     .ground_range          - geocentric ground range (Km) 
%     .group_range           - group range (Km) for each hop    
%     .phase_path            - phase path (km) for each hop
%     .initial_elev          - initial elevation (deg) at this hop   
%     .final_elev            - final elevation (deg) of this hop, IEEE NaN
%                              value is assigned if ray does not return to
%                              ground         
%     .initial_bearing       - initial ray bearing (deg) at this hop 
%     .final_bearing         - final ray bearing (deg) of this hop 
%     .deviative_absorption  - ionospheric deviative absorption (dB)    
%     .TEC_path              - total electron content along ray path (number
%                              of electrons in 1m^2 cross-section tube)  
%     .doppler_shift         - Doppler shift (Hz)
%     .geometric_path_length - Geometrical distance travelled by ray (km)
%     .frequency             - carrier frequency of the ray (MHz)
%     .nhops_attempted       - number of hops actually attempted
%     .ray_label             - label for each hop attempted which indicates
%                              what the ray has done. 
%          = 1    for ray reaching ground                           
%           -1    for field aligned backscatter - ray reflected with appropriate
%                 scattering loss, raytracing terminated
%                 *** NB FAI BACKSCATTER NOT YET IMPLEMENTED ***
%           -2    ray has penetrated the ionosphere - raytracing terminated 
%           -3    ray has exited ionospheric grid - raytracing terminated
%           -4    ray has exceeded the maximum allowed points along path
%                 (20000 points) raytracing terminated
%           -100  a catastrophic error occured - terminate raytracing
%
%   ray_path_data(:) - 1 X M structure containing information about each of the
%                      rays at each point along their paths.  The cartesian
%                      coordinate system (below) is defined with the x axis
%                      passing through the equator at the prime meridian, y
%                      through the equator at longitude of +90 degrees, and z
%                      through the geographic north pole.
%     .initial_elev            - initial elevation of the ray (degrees)
%     .initial_bearing         - initial bearing of the ray (degrees)
%     .frequency               - carrier frequency of the ray (MHz)
%     .lat                     - geocentric latitude (degrees)
%     .lon                     - geocentric longitude (degrees)
%     .height                  - height of ray (km)
%     .group_range             - group path (km)
%     .phase_path              - phase path (km)
%     .refractive_index        - refractive index
%     .group_refractive_index  - group refractive index
%     .wavenorm_ray_angle      - angle between wave normal and the ray 
%                                direction (degrees) 
%     .wavenorm_B_angle        - angle between wave normal and geomagnetic
%                                field (degrees)
%     .polariz_mag             - magnitude of the wave volume-polarization
%                                vector, R (see Notes)
%     .wave_Efield_tilt        - wave E-field tilt angle out of the plane of 
%                                the wave front 
%     .volume_polariz_tilt     - volume polarization vector tilt angle out of
%                                the wave-front plane. NB. If e- = 0 then the 
%                                magnitude of the polarization field is zero. In
%                                this case having a polarization field tilt does
%                                not make sense and so it is set IEEE NaN.
%     .electron_density        - electron density (e- / cm^3)
%     .geomag_x                - WGS84 x component of geomagnetic field (Tesla)
%     .geomag_y                - WGS84 y component of geomagnetic field (Tesla)
%     .geomag_z                - WGS84 z component of geomagnetic field (Tesla)
%     .geometric_distance      - geometrical distance travelled by ray (km)
%                                                                           
%   ray_state_vec(:) - 1 X M structure containing the state vector of each
%                      ray at each point along their paths. The cartesian
%                      coordinate system (below) is defined with the x axis
%                      passing through the equator at the prime meridian, y
%                      through the equator at longitude of +90 degrees, and z
%                      through the geographic north pole.
%      .pos_x            - cartesian x position of ray (m)
%      .pos_y            - cartesian y position of ray (m)
%      .pos_z            - cartesian z position of ray (m)
%      .dir_x            - cartesian x, y, and z components of vector 
%      .dir_y              with a magnitude equal to the refractive index  
%      .dir_z              pointing in the direction of the wave normal
%      .group_path       - group path (m)
%      .geometrical_path - geometrical distance travelled by rays (m)
%      .phase_path       - phase path (m)
%      .indep_var        - independant variable 
%      .ODE_step_size    - ODE solver independant variable step size
%                                                                            
% Outputs (2) - FOR BACKWARDS COMPATIBILITY WITH PHARLAP 3.7.1 and earlier
%   If only a single ray is being traced (i.e. length(elevs) = length(freqs) =
%   length(ray_bearings) = 1) and the environment variable PHARLAP_OLD_FORMAT
%   has been set to "true" then the output will be in the following format:
%
%   ray_data -  array containing the ray information accumulated for each hop
%     ray_data(1, :)  = geocentric end point of ray (degrees)
%     ray_data(2, :)  = geocentric longitude of end point of ray (degrees)
%     ray_data(3, :)  = geocentric ground range (Km) 
%     ray_data(4, :)  = group range (Km)                                        
%     ray_data(5, :)  = phase path (km)                                         
%     ray_data(6, :)  = initial elevation (deg) at this hop             
%     ray_data(7, :)  = final elevation (deg) of this hop, IEEE NaN value is 
%                       assigned if ray does not return to ground       
%     ray_data(8, :)  = initial ray bearing (deg) at this hop             
%     ray_data(9, :)  = final ray bearing (deg) of this hop                    
%     ray_data(10, :) = effective range (m)   *** NOT YET IMPLEMENTED *** 
%     ray_data(11, :) = ionospheric deviative absorption (dB)    
%     ray_data(12, :) = integrated electron density along ray path (number of   
%                       e- in 1m^2 cross-section tube)                         
%     ray_data(13, :) = Doppler shift (Hz)                                      
%     ray_data(14, :) = NULL (IEEE NaN) - reserved for Doppler spread (Hz) 
%     ray_data(15, :) = Geometrical distance travelled by ray (km)
%
%   ray_path_data - array containing information about the ray at each point 
%                   along it's path
%      ray_path_data(1, :) = geocentric latitude (degrees)
%      ray_path_data(2, :) = geocentric longitude (degrees)
%      ray_path_data(3, :) = height of ray above ellipsoid, km
%      ray_path_data(4, :) = group path (km)
%      ray_path_data(5, :) = phase path (km)
%      ray_path_data(6, :) = refractive index
%      ray_path_data(7, :) = group refractive index
%      ray_path_data(8, :) = angle between wave-norm and ray direction (degrees)
%      ray_path_data(9, :) = angle between wave-norm and geomag field (degrees)
%      ray_path_data(10,:) = magnitude of the wave polarization, R (see Notes)
%      ray_path_data(11,:) = wave E-field tilt angle out of the wave-front plane
%      ray_path_data(12,:) = volume polarization vector tilt angle out of the
%                            wave-front plane. NB. If e- = 0 then the 
%                            magnitude of the polarization field is zero. In 
%                            this case having a polarization field tilt does 
%                            not make sense and so it is set IEEE NaN.
%      ray_path_data(13,:) = electron density (e- / cm^3)
%      ray_path_data(14-16, :) = x,y,z components of geomagnetic field
%      ray_path_data(17, :)    = geometrical distance travelled by ray (km)
%                                                                           
%   nhops_done - number of hops actually done                               
%
%   ray_label  - array (of size nhops_done X 1) containing a label for each hop
%                which indicates what the ray has done.
%                = 1  for ray reaching ground                           
%                 -1  for field aligned backscatter - ray reflected with
%                       appropriate scattering loss, raytracing terminated
%                       *** FAI BACKSCATTER NOT YET IMPLEMENTED ***
%                 -2  ray has penetrated the ionosphere - raytracing terminated 
%                 -3  ray has exited ionospheric grid - raytracing terminated
%                 -4  ray has exceeded the maximum allowed points along path
%                       (20000 points) raytracing terminated
%                 -100  a catastrophic error occured - terminate raytracing
%
% Optional Outputs:
%   ray_state_vec - state vector of the ray at each point along its path.
%      ray_state_vec(1 - 3, :)  = x, y and z position of ray (m)
%      ray_state_vec(4 - 6, :)  = x, y, and z components of vector with a 
%                                 magnitude equal to the refractive index  
%                                 pointing in the direction of the wave normal
%     ray_state_vec_in(7, :)    = group path
%     ray_state_vec_in(8, :)    = geometrical path
%     ray_state_vec_in(9, :)    = phase path
%     ray_state_vec_in(10, :)   = independant variable 
%     ray_state_vec_in(11, :)   = ODE solver independant variable step size
%
% Notes:
% 1. Polarization :
%   The formulae to calculate the wave polarization are those from Davies (1990)
%   pp78 and assume that collisions are negligible, for reference see also 
%   Budden (1985), pp66-74. The coordinate system (cartesian i, j, k) is defined
%   as follows : the wave-normal is in the i-direction and the 
%   geomagnetic field direction lies in the i-j plane (thus k is orthoganal 
%   to B-field). Then the magnitude of the wave polarization, R, gives the
%   axial ratio of the ellipse in the j-k plane that is traced out by the 
%   electric field vector of the wave projected onto the j-k plane. The 
%   electric field vector of the wave actually lies in the plane tilted  
%   forward (or backward) by angle polariz_E_ang from the j axis for the  
%   O-mode (or X-mode). For an O-mode (or X-mode) ray, the semi-major axis of
%   the polarization ellipse projected in the j-k plane is in the j-direction
%   (or k-direction). See Figures 3.1 and 3.3 of Davies (1990).
%
% 2. Memory usage :
%   If large ionospheric grid sizes are used then the required memory could
%   exceed that available. In this case memory paging will occur which will
%   greatly slow down the computation speed. If the maximum allowable size
%   for the ionospheric grids is specified, then ~4.4GB of memory is required
%   (~1.8GB on 32-bit Windows) in addition to all the other memory
%   requirements (701 X 701 X 401 elements X 3 iono grids X 8 bytes). 
%
% References:
%   1. Haselgrove, C.B., and Haselgrove, J. (1960), "Twisted ray paths in the
%      ionosphere", Proc. Phys. Soc. (London), Vol. 75, 357-361.
%
%   2. Haselgrove, J. (1963), "The Hamiltonian ray path equations", J. Atmos. 
%      Terr. Phys., Vol. 25, 397-399. 
%  
%   3. Davies, K. (1990), "Ionospheric Radio", IEE Electromagnetic Waves 
%      Series 31, Peter Peregrinus, London.
%
%   4. Budden, K. G. (1985), "The propagation of radio waves", Cambridge 
%      University Press, Cambridge.
%
%   5. Bennett, J. A. (1967), "The calculation of Doppler Shifts due to a 
%      changing ionosphere", J. Atmos. Terr. Phys., 1967, Vol. 29, 887-891.
%

% This a Matlab help file only. The actual programme is a mex wrapper
% (raytrace-3d_matlab_wrapper.for) to the Fortran code (raytrace_3d.for).
%
% Modification history:
%   19-06-2010 M.A.Cervera  Initial version.
%                           
%   02-12-2015 M.A.Cervera
%      Updated for multiple ray input/ouput capabilities of PHaRLAP 4.0.0
%                           
