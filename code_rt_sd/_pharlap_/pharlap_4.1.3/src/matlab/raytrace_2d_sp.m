%
% Name :
%   raytrace_2d_sp
% 
% Purpose / description:
%   2D numerical raytrace for a multihop ray with correction for geomagnetic 
%   field. Geometric loss, deviative ionospheric absorption, backscattered
%   loss due to field aligned irregularites, ray becoming evanescent, 
%   Doppler shift and Doppler spread are considered. The ray state vector may
%   be returned if required. A user defined starting state vector may be
%   input. However, use this feature with caution, it is for experts
%   only. Spherical earth coordinate system is assumed.
%    
% Calling sequence:
%   The first time raytrace_2d is called, ionospheric grids
%   must be specified. For subsequent calls if these grids are not passed in 
%   (see calling sequence examples 3 and 4) then the grids from the
%   previous call (which are in memory) will be used. This can speed up the 
%   subsequent calls to raytrace for large grids. 
%
%   1. [ray_data, ray_path_data, ray_state_vec] = ...
%         raytrace_2d_sp(elev, bearing, freq, nhops, tol, radius_earth, ...
%              irregs_flag, iono_en_grid, iono_en_grid_5, collision_freq, ...
%              start_height, height_inc, range_inc, irreg);
%
%   2. [ray_data, ray_path_data, ray_state_vec] = ...
%         raytrace_2d_sp(elev, bearing, freq, nhops, tol, radius_earth, ...
%              irregs_flag, iono_en_grid, iono_en_grid_5, collision_freq, ...
%              start_height, height_inc, range_inc, irreg, ray_state_vec_in);
%
%   3. [ray_data, ray_path_data, ray_state_vec] = ...
%         raytrace_2d_sp(elev, bearing, freq, nhops, tol, radius_earth, ...
%              irregs_flag);
%
%   4. [ray_data, ray_path_data, ray_state_vec] = ...
%         raytrace_2d_sp(elev, bearing, freq, nhops, tol, radius_earth, ...
%              irregs_flag, ray_state_vec_in);
%
% Calling sequence (2) : THESE ARE FOR BACKWARDS COMPATIBILITY WITH PHARLAP 
%                        VERSION 3.7.1 AND EARLIER
%   If only a single ray is being traced (i.e. length(elevs) = 1) and the 
%   environment variable PHARLAP_OLD_FORMAT_2D has been set to "true" then the 
%   calling sequence will be:
%
%   5. [ray_data, ray_path_data, nhops_attempted, ray_label, ray_state_vec] = ...
%         raytrace_2d_sp(elev, bearing, freq, nhops, tol, radius_earth, ...
%              irregs_flag, iono_en_grid, iono_en_grid_5, collision_freq, ...
%              start_height, height_inc, range_inc, irreg);
%
%   6. [ray_data, ray_path_data, nhops_attempted, ray_label, ray_state_vec] = ...
%         raytrace_2d_sp(elev, bearing, freq, nhops, tol, radius_earth, ...
%              irregs_flag, iono_en_grid, iono_en_grid_5, collision_freq, ...
%              start_height, height_inc, range_inc, irreg, ray_state_vec_in);
%
%   7. [ray_data, ray_path_data, nhops_attempted, ray_label, ray_state_vec] = ...
%         raytrace_2d_sp(elev, bearing, freq, nhops, tol, radius_earth, ...
%              irregs_flag);
%
%   8. [ray_data, ray_path_data, nhops_attempted, ray_label, ray_state_vec] = ...
%         raytrace_2d_sp(elev, bearing, freq, nhops, tol, radius_earth, ...
%              irregs_flag, ray_state_vec_in);
%
% Inputs (scalar unless indicated):
%   elevs          - 1 X M (where M is the number of rays) array of initial
%                    elevation of rays (deg) 
%   bearing        - bearing (degrees from North) of ray
%   freq           - wave frequency (MHz)    
%   nhops          - number of hops to be performed
%   tol            - a 1 or 3 element vector controlling ODE solver precision.
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
%   radius_earth   - radius of the Earth (Km) to use for spherical Earth 
%                      coordinate system
%   irreg_flag     - flag indicating if irregularities are to be turned on (= 1)
%                    or not (= 0). If the irregularities are turned off then 
%                    (1) field aligned backscatter turned off and (2) Doppler 
%                    spread is set to 0.  
%   iono_en_grid   - 2d grid (height vs ground range) of ionospheric 
%                       electron density (electrons / cm^3). Maximum grid
%                       size is 2001 X 2001 elements.
%   iono_en_grid_5 - 2d grid (height vs ground range) of ionospheric 
%                       electron density (electrons / cm^3) 5 minutes later.
%                       Maximum grid size is 2001 X 2001 elements.
%   collision_freq - 2d grid (height vs ground range) of collision freqs (Hz).
%                       Maximum grid size is 2001 X 2001 elements.
%   start_height   - start height of iono_en_grid and bfield grids (km)
%   height_inc     - height step of iono_en_grid and bfield grids (km)
%   range_inc      - range step of iono_en_grid, bfield, iono_parms, dec,
%                      dip, and irreg_strngth arrays (km)
%   irreg         -  4 x num_ranges array of irregularity parameters as a
%                       function of ground range. Maximum grid size is 
%                       4 X 2001 elements. 
%     irreg(1, :) =  irregularity strength and is the ratio of irregular
%                      electron density to the background value - can be
%                      ignored (set to 0) if irreg_flag = 0
%     irreg(2, :) =  magnetic dip angle at the phase screen height of 
%                      irregularities (typically 300km) (degrees) - can be 
%                      ignored (set to 0) if irreg_flag = 0
%     irreg(3, :) =  magnetic declination at the phase screen height of 
%                      irregularities (typically 300km) (degrees) - can be 
%                      ignored (set to 0) if irreg_flag = 0
%     irreg(4, :) =  square of frequency spread (Hz^2) per unit path length (Km)
%                      at a carrier frequency of 1MHz scaled by the electron
%                      density (cm^-3) - required for Doppler spread calculation
%                      - can be ignored (set to 0) if irreg_flag = 0
%
% Optional Inputs:
%   ray_state_vec_in(:) - 1 X M structure containing the state vector of each
%            ray (total of M rays) at its starting point. Each field (total of
%            9 fields) is an element of the state vector. If the first field
%            for a given ray is -1 then default starting values are used for
%            that ray. NB: non-default values are for EXPERTS ONLY - you really
%            want to know what you are doing!
%     .r                     - Distance of ray to centre of Earth (km)
%     .Q                     - This is Q (eqn 4. of Coleman JASTP, 59, pp2090). 
%                              At ground level its value is sin(ray elevation)
%     .theta                 - Angle subtended by ray at the centre of the Earth
%                              (radians)
%     .delta_r               - delta_r (see eqn 7 of Coleman RS, 33, pp1188).
%                              Required to calculate focussing gain/loss
%     .delta_Q               - delta_Q (see eqn 8 of Coleman RS, 33, pp1188). 
%                              Required to calulate focussing gain/loss
%     .deviative_absorption  - Ionospheric deviative absorption  (dB)     
%     .phase_path            - Phase path (km) 
%     .group_path            - Group path (km) - independant var for ODE solver
%     .group_path_step_size  - Group path step size (km) for RKF ode solver
%
%
% Outputs:
%   ray_data(:) -  1 X M structure containing the information for each of the
%                  rays (total of M rays). Each field is a 1 X N array
%                  containing information for each hop. The fields of this
%                  structure are:
%     .ground_range          - geodetic (WGS84) ground range (Km)
%     .group_range           - group range (Km)    
%     .phase_path            - phase path (km)
%     .geometric_path_length - physical distance along ray path (Km)
%     .initial_elev          - initial elevation (deg) of this hop
%     .final_elev            - final elevation (deg) of this hop
%     .apogee                - maximum altitude of ray (Km) 
%     .gnd_rng_to_apogee     - ground range to max height (km)  
%     .plasma_freq_at_apogee - plasma frequency at the ray's apogee (MHz)
%     .virtual_height        - virtual height (km)
%     .effective_range       - effective range (m)
%     .deviative_absorption  - ionospheric deviative absorption(dB)
%     .TEC_path              - integrated electron density along ray path
%                              (number of electrons in 1m^2 cross-section tube)
%     .Doppler_shift         - Doppler shift (Hz)
%     .Doppler_spread        - Doppler spread (Hz)              
%     .FAI_backscatter_loss  - backscattered loss (dB) for the last hop due to 
%                              field aligned irregularites.If there are no FAIs 
%                              then this is set to 0 for all hops.
%     .nhops_attempted       - number of hops actually attempted
%     .ray_label             - label for each hop attempted which indicates
%                              what the ray has done. 
%           = 1  for ray reaching ground                           
%             0  for ray becoming evanescent, raytracing terminated
%            -1  for field aligned backscatter - ray reflected with
%                appropriate scattering loss, raytracing terminated
%            -2  ray has penetrated the ionosphere - raytracing terminated 
%            -3  ray has exceeded max. ground range - raytracing terminated
%            -4  ray angular coordinate has become negative (bad - should 
%                never happen) - raytracing terminated
%            -5  ray has exceeded the maximum allowed points along path
%                (20000 points) - raytracing terminated
%            -6  ray is near antipodal point, the WGS84 coordinate
%                conversion routines are unreliable - terminate
%                raytracing 
%          -100  a catastrophic error occured - terminate raytracing
%
%   ray_path_data(:) - 1 X M structure containing information about each of the
%                      rays at each point along their paths 
%     .ground_range       - geodetic (WGS84) ground range from origin to 
%                           point on ground directly below ray, km
%     .height             - height of ray above WGS84 ellipsoid, km
%     .group_range        - group range, km
%     .phase_path         - phase path, km
%     .geometric_distance - physical distance along ray path, km
%     .electron_density   - electron density, 1/cm^3
%     .refractive_index   - refractive index
%
%   ray_state_vec(:) - 1 X M structure containing the state vector of each
%                      ray at each point along their paths.
%     .r                    - Distance of ray to centre of Earth (km)
%     .Q                    - This is Q (eqn 4. of Coleman JASTP, 59, pp2090)
%                             at ground level its value is sin(ray elevation)
%     .theta                - Angle subtended by ray at the centre of the Earth
%                             (radians)
%     .delta_r              - delta_r (km) (see eqn 7 of Coleman RS, 33, pp1188)
%                             Required to calculate focussing gain/loss
%     .delta_Q              - delta_Q (see eqn 8 of Coleman RS, 33, pp1188).
%                             Required to calulate focussing gain/loss
%     .deviative_absorption - Ionospheric deviative absorption  (dB)     
%     .phase_path           - Phase path (km)  
%     .group_path           - independant var for RKF ode solver
%     .group_step_size      - Group path step size (km) for  ode solver
%
% Outputs (2) - FOR BACKWARDS COMPATIBILITY WITH PHARLAP 3.7.1 and earlier
%   If only a single ray is being traced (i.e. length(elevs) = length(freqs) =
%   length(ray_bearings) = 1) and the environment variable PHARLAP_OLD_FORMAT
%   has been set to "true" then the output will be in the following format:
%   ray_data   - array (of size 19 X nhops_attempted) containing the ray
%                 information accumulated for each hop :
%     ray_data(1, :)  = NULL (IEEE NaN)
%     ray_data(2, :)  = NULL (IEEE NaN)
%     ray_data(3, :)  = geodetic (WGS84) ground range (Km)
%     ray_data(4, :)  = group range (Km)    
%     ray_data(5, :)  = maximum height (Km) 
%     ray_data(6, :)  = ground range to max height (km)  
%     ray_data(7, :)  = initial elevation (deg) of this hop
%     ray_data(8, :)  = final elevation (deg) of this hop
%     ray_data(9, :)  = Doppler spread (Hz)              
%     ray_data(10, :) = Doppler shift (Hz)
%     ray_data(11, :) = physical distance along ray path (Km)
%     ray_data(12, :) = effective range (m)
%     ray_data(13, :) = ionospheric deviative absorption(dB)
%     ray_data(14, :) = plasma frequency at the ray's apogee (MHz)
%     ray_data(15, :) = backscattered loss (dB) for the last hop due to field 
%                         aligned irregularites.If there are no FAIs then
%                         this is set to 0 for all hops.
%     ray_data(16, :) = virtual height (km)
%     ray_data(17, :) = phase path (km)
%     ray_data(18, :) = integrated electron density along ray path  (number of
%                         electrons in 1m^2 cross-section tube)
%     ray_data(19, :) = NULL (IEEE NaN) - reserved for future use 
%
%   ray_path_data - array containing information about the rays's path
%     ray_path_data(1, *) = geodetic (WGS84) ground range from origin to 
%                            point on ground directly below ray, km
%     ray_path_data(2, *) = height of ray above WGS84 ellipsoid, km
%     ray_path_data(3, *) = group range, km
%     ray_path_data(4, *) = phase path, km
%     ray_path_data(5, *) = physical distance along ray path, km
%     ray_path_data(6, *) = electron density, 1/cm^3
%     ray_path_data(7, *) = refractive index
%
%   nhops_attempted - number of hops actually attempted
%
%   label      - array (of size nhops_attempted X 1) containing a label for each 
%                hop which indicates what the ray has done.
%                = 1  for ray reaching ground                           
%                  0  for ray becoming evanescent, raytracing terminated
%                 -1  for field aligned backscatter - ray reflected with
%                       appropriate scattering loss, raytracing terminated
%                 -2  ray has penetrated the ionosphere - raytracing terminated 
%                 -3  ray has exceeded max. ground range - raytracing terminated
%                 -4  ray angular coordinate has become negative (bad - should 
%                     never happen) - raytracing terminated
%                 -5  ray has exceeded the maximum allowed points along path
%                       (20000 points) - raytracing terminated
%                 -100  a catastrophic error occured - terminate raytracing
%
% Optional Outputs:
%   ray_state_vec -  state vector of the ray at each point along the ray path
%     ray_state_vec(1, *) = Distance of ray to centre of Earth (km)
%     ray_state_vec(2, *) = This is Q (eqn 4. of Coleman JASTP, 59, pp2090)    
%                             at ground level its value is sin(ray elevation)
%     ray_state_vec(3, *) = Angle subtended by ray at the centre of the Earth
%                             (radians)
%     ray_state_vec(4, *) = delta_r (km) (see eqn 7 of Coleman RS, 33, pp1188).
%                             Required to calculate focussing gain/loss
%     ray_state_vec(5, *) = delta_Q (ess eqn 8 of Coleman RS, 33, pp1188).     
%                             Required to calulate focussing gain/loss
%     ray_state_vec(6, *) = Ionospheric deviative absorption  (dB)     
%     ray_state_vec(7, *) = Phase path (km)                                     
%     ray_state_vec(8, *) = Group path (km) - independant var for RKF ode solver
%     ray_state_vec(9, *) = Group path step size (km) for RKF ode solver
%


% This a Matlab help file only. The actual programme is a mex wrapper
% (raytrace-2d-sp_matlab_wrapper.for) to the Fortran code (raytrace_2d.for).