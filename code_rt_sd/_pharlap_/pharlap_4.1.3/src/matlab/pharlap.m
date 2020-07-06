%
% PHaRLAP toolbox documentation
%
% Name :
%   PHaRLAP - Provision of High-frequency Raytracing LAboratory for Propagation
%             studies
%
% Purpose :
%   Provides 2D and 3D numerical raytracing engines and supporting routines 
%   (geomagnetic, ionospheric, etc.) to enable modeling of radio wave
%   propogation through the ionosphere. Geometric focussing, ionospheric
%   absorption (D-region via Bradley-George, above D-region and deviative),
%   ground forward scatter and backscatter losses, backscattered loss due to
%   field aligned irregularites, ray becoming evanescent, O-X splitting,
%   polarization, Doppler shift and Doppler spread are considered. Raytracing 
%   through Es patches is not considered (however see ray_test2.m as a
%   potential method for modelling Es). The ray trajectories and the full
%   state-vector of the ray at each  point along the ray trajectory is
%   available. User modified state-vectors may be input to the raytracing
%   routines for advanced ray studies. The coordinate system used is WGS84. 
%
% Supplied routines (type help <routine name> for further help) :
%   
%   Core routines
%   -------------
%   raytrace_2d         - 2D numerical raytracing engine for a multihop ray 
%                         (assumes WGS84 coordinate system)
%   raytrace_2d_sp      - 2D numerical raytracing engine for a multihop ray 
%                         (assumes spherical Earth coordinate system).
%   raytrace_2d_gm      - uses raytrace_2d and applies O-X correction factor
%                         calculated from gm_freq_offset so that approximations
%                         for O and X mode 2D rays may be obtained
%   raytrace_3d         - full 3D magneto-iono numerical raytrace engine for
%                         O-mode, X-mode and "no-field" rays (assumes WGS84
%                         coordinate system)
%   raytrace_3d_sp      - full 3D magneto-iono numerical raytrace engine for
%                         O-mode, X-mode and "no-field" rays (assumes
%                         spherical Earth coordinate system with user
%                         specified radius)
%   abso_bg             - D-region absorption via Bradley-George
%   chapman             - calculates ionospheric plasma profiles based on
%                         Chapman layers
%   dop_spread_eq       - simple model of Doppler spread imposed on ray 
%                         traversing the equatorial region
%   ground_bs_loss      - power loss of the radio-waves back-scattered from
%                         the ground
%   ground_fs_loss      - calculates forward ground scattering losses
%   gen_iono_grid_2d    - generates ionospheric parameters array, ionospheric 
%                         plasma density grid, and irregularity strength in
%                         format required by the 2D raytracing engine - uses
%                         iri2007
%   gen_iono_grid_3d    - generates ionospheric parameters array, ionospheric 
%                         plasma density grid, geomagnetic field grid, and 
%                         irregularity strength in format required by the 3D 
%                         raytracing engine - uses iri2007 and igrf2007  
%   gm_freq_offset      - calculates the approximate geomagnetic O-X mode
%                         frequency split (MHz) for a specified propagation 
%                         path.
%   igrf2007            - International Geomagnetic Reference Field (distributed
%                         with IRI2007).
%   igrf2011            - International Geomagnetic Reference Field (distributed
%                         with IRI2012).
%   igrf2016            - International Geomagnetic Reference Field (distributed
%                         with IRI2016).
%   iri2007             - International Reference Ionosphere (2007)
%   iri2012             - International Reference Ionosphere (2012) 
%   iri2016             - International Reference Ionosphere (2016) 
%   iri2012_firi_interp - calls International Reference Ionosphere (2012)
%                         with FIRI rocketsonde-based lower ionosphere and
%                         performs interpolation/smooth to remove discontinuity
%   iri2016_firi_interp - calls International Reference Ionosphere (2016)
%                         with FIRI rocketsonde-based lower ionosphere and
%                         performs interpolation/smooth to remove discontinuity
%   irreg_strength      - simple model of irregularity strength in the
%                         equatorial and auroral regions
%   land_sea            - returns land or sea for given location on Earth
%   nrlmsise00          - NRLMSISe-00 model atmos. (distributed with IRI2012)
%   noise_ccir          - CCIR environmental noise model
%   plot_ray_iono_slice - plots ionospheric slice in an arc (to preserve
%                         Earth geometry) and overplots rays
%
%   Ancilliary routines
%   -------------------
%   coning              - calulates correction to azimuth of ray due to the  
%                         cone effect of linear arrays 
%   deriv               - calculates derivative via 3-point, Lagrangian 
%                         interpolation
%   earth_radius_wgs84  - returns the WGS84 radius of the Earth at input 
%                         geodetic latitude
%   ENU2xyz             - convert the location of a point specified in an East, 
%                         North, Up frame at a local origin on the Earth to
%                         cartesian coordinates (x, y, z) 
%   julday              - calculates the Julian Day Number 
%   latlon2raz          - converts spherical Earth and longitude to range and 
%                         azimuth with respect to an origin for various geoids
%   raz2latlon          - converts range and azimuth from origin to spherical 
%                         Earth latitude and longitude for various geoids
%   solar_za            - returns the solar zenith angle
%   wgs842gc_lat        - convert WGS84 geodetic latitude to geocentric
%                         latitude  
%   wgs84_llh2xyz       - converts WGS64 geodetic lat., long. and height to 
%                         Earth centred x, y, z coordinates
%   wgs84_xyz2llh       - converts Earth centred x, y, z coordinates to WGS84
%                         geodetic lat., long. and height
%   xyz2ENU             - convert the location of a point specified in an
%                         cartesian coordinate (x, y, z) frame at a local
%                         origin on the Earth to East, North, Up coordinates.
%
%   Examples
%   --------
%   ray_test1           - simple 2D NRT example showing fan of rays using
%                         WGS84 coordinates
%   ray_test2           - 2D NRT example showing ray state vector modification
%   ray_test3           - simple 2D NRT showing fan of rays using spherical
%                         Earth corrdinates
%   ray_test4           - simple 2D NRT showing rays with different carrier
%                         frequencies
%   ray_test_3d         - 3D NRT (WGS84 coordinates) example showing O, X and 
%                         no-field rays
%   ray_test_3d_sp      - 3D NRT (spherical Earth coordinates) example
%                         showing O, X and no-field rays
%   ray_test_3d_pol     - 3D NRT for single ray with polarization calculations 
%   ray_test_3d_iono_tilt - 3D NRT thorough a "tilted" ionosphere
%   ois_synth           - example of synthetic single hop OIS ionogram
%                         generation using 2D numerical raytracing
%   ois_synth_mh        - multi-hop OIS ionogram synthesis with GUI
%   bss_synth           - example of multi-hop back-scatter ionogram synthesis
%   NRT_comparision     - comparison of a fan of rays using 2D and 3D NRT
%
% Notes:
%   1. The routines irreg_strength and dop_spread_eq are required for 
%      propagation studies where Doppler spread is desired. Currently
%      they are based on very simple models. Treat the these results with
%      caution. See M.A. Cervera for further details. 
%
% Author:
%   M.A. Cervera  16/06/2006
%   Last update to this file:  07/07/2016 (M.A. Cervera)
%
% Modification History:
%   See release notes.
%



%  This file constitutes documentation only.
