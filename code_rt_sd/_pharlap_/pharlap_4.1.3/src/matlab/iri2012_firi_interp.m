%
% Name :
%   iri2012_firi_interp.m
%
% Purpose :
%   Matlab wrapper to the IRI-2012 fortan based empirical model ionosphere,
%   using the FIRI rocketsonde-based model for the region below 120 km, and
%   interpolating between the FIRI model and the default IRI F layer model.
%   For (1) a given user specified R12 index (i.e. ionospheric condition), or 
%   (2) a given epoch, or (3) user specified ionospheric layers, it returns
%   modeled electron density, electron temperature, ion temperature, and ion
%   composition (O+, H+, He+, NO+, O+2) height profiles, ion drift and F1
%   probability along with various other parameters including foE, hmE, foF1,
%   hmF1, foF2, hmF2. Storms are not modelled for case (1) and may be
%   switched on or off for case (2). NB: SEE NOTES SECTION BELOW FOR
%   IMPORTANT INFORMATION.  
%
% Calling sequence :
%    1. [iono, iono_extra] = iri2012_firi_interp(lat, lon, R12, UT, ...
%                                 ht_start, ht_step, num_hts);
%
%    2. [iono, iono_extra] = iri2012_firi_interp(lat, lon, R12, UT, ...
%                                 ht_start, ht_step num_hts, 
%                                 'iono_layer_parms', iono_layer_parms);
%
%    3. [iono, iono_extra] = iri2012_firi_interp(lat, lon, R12, UT, ...
%                                 ht_start, ht_step num_hts, ...
%                                 'B0B1_model', B0B1_model, ...
%                                 'iono_layer_parms', iono_layer_parms);
%
% Inputs :
%   lat      - geographic latitude of point (degrees)  
%   lon      - geographic longitude of point (degrees)
%   R12      - depends on value as follows:
%      R12 = 1 - 200 :  IRI2012 is called with R12 input as the user specified
%                       yearly smoothed monthly median sunspot number. The
%                       storm model is turned off
%      R12 = -1      :  IRI2012 is called with ionospheric conditions read
%                       from file (ig_rz.dat) based on input epoch (UT) and
%                       may be historical or projected conditions (dependent
%                       on the epoch). See notes below.
%      R12 = -2      :  as for R12 = -1 but storm model is turned off 
%
%   UT       - 5x1 array containing UTC date and time - year, month, day, 
%              hour, minute
%   ht_start - start height for ionospheric profile generation (km)
%   ht_step  - height step size for ionospheric profile generation (km)
%   num_hts  - number of heights for ionospheric profile generation
%              must be >= 2 and <= 1000
%
% Optional Inputs:
%   B0B1_model - specifies the model to use for the bottomside ionospheric
%                profile. The model defines how IRI2012 determines the B0
%                (thickness) and B1 (shape) parameters. If B0B1_model is not 
%                specified then the default is Bil-2000 (B0B1_model = 2)
%       B0B1_model = 1 : ABT-2009  (Adv. Space Res., Vol 43, 1825-1834, 2009)
%       B0B1_model = 2 : Bil-2000  (Adv. Space Res., Vol 25, 89-96, 2000)
%       B0B1_model = 3 : Gul-1987  (Adv. Space Res., Vol 7, 39-48, 1987)
%
%   iono_layer_parms - 6 element array containing user specified ionospheric 
%                      layer parameters. Specification of this optional input
%                      overrides the model parameters in IRI. Critical
%                      frequencies must be in the range 0.1MHz to 100MHz with
%                      foE < foF1 < foF2. Layer heights must be in the range
%                      50km to 1000km with hmE < hmF1 < hmF2. Parameters with
%                      a value set to -1 indicate to IRI2012 to use its model
%                      for that case. 
%       iono_layer_parms(1) = foF2
%       iono_layer_parms(2) = hmF2
%       iono_layer_parms(3) = foF1
%       iono_layer_parms(4) = hmF1
%       iono_layer_parms(5) = foE
%       iono_layer_parms(6) = hmE
%   
% Outputs :
%   iono  -  Array of 11 output parameters x num_hts heights. NB, ionosperic
%            profiles are not computed below 65 km or above 2000 km. Values
%            of -1 are returned if these heights are requested and also for
%            those heights where a valid value is unable to be calculated. If
%            the optional inputs  ht_start and ht_step are not supplied then
%            the profiles are not calculated.
%       iono(1, :) = electron number density (m^-3)
%       iono(2, :) = neutral temperature (K)
%       iono(3, :) = ion temperature (K)
%       iono(4, :) = electron temperature (K)
%       iono(5, :) = O+ ion density (m^-3)
%       iono(6, :) = H+ ion density (m^-3)
%       iono(7, :) = He+ ion density (m^-3)
%       iono(8, :) = O2+ ion density (m^-3)
%       iono(9, :) = NO+ ion density (m^-3)
%       iono(10, :) = cluster ions density (m^-3)
%       iono(11, :) = N+ ion density (m^-3)
%
%   iono_extra - Array of extra output parameters. NB, Unused array elements  
%                and parameters which are not able to be calculated are flagged 
%                by a value of -1.
%      iono_extra(1) = NmF2 (m^-3)         iono_extra(2) = HmF2 (km)
%      iono_extra(3) = NmF1 (m^-3)         iono_extra(4) = HmF1 (km)
%      iono_extra(5) = NmE (m^-3)          iono_extra(6) = HmE (km)
%      iono_extra(7) = NmD (m^-3)          iono_extra(8) = HmD (km)
%      iono_extra(9) = HHALF (km)          iono_extra(10) = B0 (km)
%      iono_extra(11) = VALLEY-BASE (m^-3) iono_extra(12) = VALLEY-TOP (km)
%      iono_extra(13) = Te-PEAK (K)        iono_extra(14) = Te-PEAK HEIGHT (km)
%      iono_extra(15) = Te-MOD(300KM) (K)  iono_extra(16) = Te-MOD(400KM) (K)
%      iono_extra(17) = Te-MOD(600KM) (K)  iono_extra(18) = Te-MOD(1400KM) (K)
%      iono_extra(19) = Te-MOD(3000KM) (K) iono_extra(20) = Te(120KM)=TN=Ti (K)
%      iono_extra(21) = Ti-MOD(430KM) (K)  iono_extra(22) = X (km), where Te=Ti
%      iono_extra(23) = sol zen. ang (deg) iono_extra(24) = sun declin. (deg)
%      iono_extra(25) = DIP (deg)          iono_extra(26) = dip latitude (deg)
%      iono_extra(27) = modified dip lat.  iono_extra(28) = DELA
%      iono_extra(29) = sunrise (hours)    iono_extra(30) = sunset (hours)
%      iono_extra(31) = ISEASON (1=spring) iono_extra(32) = NSEASON (northern)
%      iono_extra(33) = Rz12               iono_extra(34) = Covington Index
%      iono_extra(35) = B1                 iono_extra(36) = M(3000)F2
%      iono_extra(37) = Unused             iono_extra(38) = Unused
%      iono_extra(39) = gind (IG12)        iono_extra(40) = F1 probability (old)
%      iono_extra(41) = F10.7 daily        iono_extra(42) = c1 (F1 shape)
%      iono_extra(43) = daynr              iono_extra(44) = equatorial vertical 
%      iono_extra(45) = foF2_storm/foF2_quiet               ion drift in m/s
%      iono_extra(46) = F10.7_81           iono_extra(47) = foE_storm/foE_quiet 
%      iono_extra(48) = spread-F probability          
%      iono_extra(49) = Geomag. latitude   iono_extra(50) = Geomag. longitude  
%      iono_extra(51) = Ap at current time iono_extra(52) = daily Ap
%      iono_extra(53) = invdip/degree      iono_extra(54) = MLT
%      iono_extra(55) = CGM-latitude       iono_extra(56) = CGM-longitude
%      iono_extra(57) = CGM-lati(MLT=0)    iono_extra(58) = CGM-lati for MLT=1
%      iono_extra(59) = CGM-lati(MLT=2)    iono_extra(60) = CGM-lati for MLT=3
%      iono_extra(61) = CGM-lati(MLT=4)    iono_extra(62) = CGM-lati for MLT=5
%      iono_extra(63) = CGM-lati(MLT=6)    iono_extra(64) = CGM-lati for MLT=7
%      iono_extra(65) = CGM-lati(MLT=8)    iono_extra(66) = CGM-lati for MLT=9
%      iono_extra(67) = CGM-lati(MLT=10)   iono_extra(68) = CGM-lati for MLT=11
%      iono_extra(69) = CGM-lati(MLT=12)   iono_extra(70) = CGM-lati for MLT=13
%      iono_extra(71) = CGM-lati(MLT=14)   iono_extra(72) = CGM-lati for MLT=15
%      iono_extra(73) = CGM-lati(MLT=16)   iono_extra(74) = CGM-lati for MLT=17
%      iono_extra(75) = CGM-lati(MLT=18)   iono_extra(76) = CGM-lati for MLT=19
%      iono_extra(77) = CGM-lati(MLT=20)   iono_extra(78) = CGM-lati for MLT=21
%      iono_extra(79) = CGM-lati(MLT=22)   iono_extra(80) = CGM-lati for MLT=23
%      iono_extra(81) = CGM-MLT            iono_extra(82) = CGM-lati for CGM-MLT
%      iono_extra(83) = Kp at current time iono_extra(84) = magnetic declination
%      iono_extra(85) = L-value            iono_extra(86) = dipole moment 
%      iono_extra(87 - 100) = Unused
%
%  Notes :
%   1. Notes for IRI2012 called using specified input ionospheric conditions:
%   1.1 If the ionospheric conditions are controlled by the matlab input R12
%       index then the input year (in the UT array) has a very small effect
%       on  the solar conditions. For example 
%          [iono iono_extra] = iri2012(-25, 135, 70, [2000 1 1 3 0])
%       returns NmF2 = 1.0252e+12 electrons/m^3, whereas
%          [iono iono_extra] = iri2012(-25, 135, 70, [2001 1 1 3 0])
%       returns NmF2 = 1.0260e+12
%
%   1.2 User defined IG12, F10.7 and F10.7_81 (3 solar rotation average of
%       F10.7) required by IRI2012 are derived from R12 using the following
%       empirical formulas : 
%            F107    = 63.7 + 0.728*R12 + 0.00089*R12^2
%            F107_81 = F107
%            IG12 = -12.349154 + R12 * (1.4683266 - R12 * 2.67690893e-03);
%       These derived values for IG12, F10.7 and F10.7_81 are input into 
%       IRI-2012
%
%   2. Notes for IRI2012 called using specified input epoch:
%   2.1 IRI2012 uses solar indices tabled in ig_rz.dat (which is supplied with
%       IRI2012). This file contains R12 and IG12 from Jan 1958 to Dec 2018. If 
%       the input UT is outside this range then an error is returned. R12 is
%       used to model the height of the F2 layer and IG12 its strength.
%
%   2.2 The computation of the yearly-running mean for month M requires the
%       indices for the six months preceeding M and the six months following 
%       M (month: M-6, ..., M+6). To calculate the current running mean one 
%       therefore requires predictions of the index for the next six months. 
%       Starting from six months before the UPDATE DATE (listed at the top of 
%       the file ig_rz.dat and below at point 3.) and onward the indices are 
%       therefore based on indices predictions.
%
%   2.3 ig_rz.dat updated 20-Feb-2015
%
%   2.4 The solar activity parameter F10.7 (daily) and magnetic activity Ap
%       index (3-hourly) used by IRI2012 are tabled in apf107.dat from 1 Jan 
%       1958 to 31 Dec 2014. If UT is outside this range then the storm model
%       (which relies on  Ap) is switched off and a monthly median F10.7
%       (calculated from R12 using the empirical formula F107 = 63.7 +
%       0.728*R12 + 0.00089*R12^2) is used in  place of the daily
%       F10.7. 
%
%   3. This mex file drives IRI-2012 with the following options
%      3.1  Ne computed
%      3.2  Te, Ti computed
%      3.3  Ne & Ni computed
%      3.4  B0,B1 - Bil-2000 (unless over-ridden by user)       
%      3.5  foF2  - URSI 
%      3.6  Ni    - RBV-2010 & TTS-2005
%      3.7  Ne    - f10.7 unlimited
%      3.8  foF2 from model (unless over-ridden and input by user)
%      3.9  hmF2 from model (unless over-ridden and input by user)
%      3.10 Te    - Standard
%      3.11 Ne    - Standard Profile
%      3.12 Messages to unit 6 (but see 3.34 below)
%      3.13 foF1 from model (unless over-ridden and input by user)
%      3.14 hmF1 from model (unless over-ridden and input by user)
%      3.15 foE  from model (unless over-ridden and input by user)
%      3.16 hmE  from model (unless over-ridden and input by user)
%      3.17 Rz12 from file (unless over-ridden and input by user, i.e. R12 > 0)
%      3.18 IGRF magnetic field model
%      3.19 F1 probability model
%      3.20 standard F1
%      3.21 ion drift computed
%      3.22 ion densities in m-3
%      3.23 Te_topside - TBT-2012
%      3.24 D-Region model - IRI-1990 (unless over-ridden by user)
%      3.25 F107D from APF107.DAT (unless over-ridden by user, i.e. R12 > 0)
%      3.26 foF2 no storm updating (unless over-ridden by user)
%      3.27 IG12 from file (unless over-ridden and input by user, i.e. R12 > 0)
%      3.28 spread-F probability - not computed
%      3.29 false i.e. topside defined by 3.30 below
%      3.30 NeQuick topside model 
%      3.31 not used as 3.4 is set to Bil-2000 (unless over-ridden by user)
%      3.32 F10.7_81 from file  (unless over-ridden by user, i.e. R12 > 0)
%      3.33 Auroral boundary model is off
%      3.34 Messages off
%      3.35 no foE storm updating
%      3.36 hmF2 without foF2-storm                
%      3.37 topside without foF2-storm 
%      3.38 turn WRITEs off in IRIFLIP
%
% Further information :
%   http://iri.gsfc.nasa.gov/
%   http://irimodel.org/    
%
% Modification History:
%   V1.0  L. H. Pederick 03/12/2012
%
%   V1.1  L. H. Pederick 09/09/2013
%      Fixed bug where interpolation zone was positioned incorrectly
%      (top_transition < bottom_transition)
%
%   V1.2  M. A. Cervera  16/12/2013
%      Now accepts user input of foF2, hmF2, foF1, hmF1, foE and hmE
%      if required
%
%   V1.3  M. A. Cervera  30/06/2015
%      Modified to be consistent with new version of iri2012
%

function [iono, iono_extra] = iri2012_firi_interp(site_lat, site_lon, R12, ...
      UT, start_height, height_step, num_heights, varargin)
  
  if (nargin < 7)
    error('iri2012_firi_interp:argChk', ...
	  'Wrong number of input arguments: at leaset 7 inputs required')
    return
  end

  iri2012_arg_str = ['site_lat, site_lon, R12, UT, start_height, ' ...
		  'height_step, num_heights, B0B1_model, D_model'];
  
  D_model = 2;     % selects FIRI
  B0B1_model = 2;  % the default B0B1_model used for IRI
  
  for ii = 1:nargin-7
    % if B0B1_model has been input then over-ride the default
    if strcmp(varargin{ii}, 'B0B1_model')  
      if nargin-7 < ii+1
	error('B0B1_model not set')
      end
      B0B1_model = varargin{ii+1};
      if length(B0B1_model) > 1
        error('invalid value for B0B1_model')
      end
      if ~isreal(B0B1_model) | ischar(B0B1_model) | isnan(B0B1_model) | ...
	 ~isfinite(B0B1_model)
        error('invalid value for B0B1_model')
      end     
    end
    
    % if iono_layer_parms has been input then add it to argument list for iri
    if strcmp(varargin{ii}, 'iono_layer_parms')
      if nargin-7 < ii+1
	error('iono_layer_parms not set')
      end
      if length(varargin{ii+1}) ~=6
	error('invalid input for iono_layer_parms')
      end
      iono_parms = varargin{ii+1}(1:6);
      if ~isreal(iono_parms) | ischar(iono_parms) | isnan(iono_parms) | ...
	 any(~isfinite(iono_parms))
        error('invalid value for iono_parms')
      end           
      iri2012_arg_str = [iri2012_arg_str ', iono_parms'];
    end
    
  end   

  iricall_str = ['[iono iono_extra] = iri2012(' iri2012_arg_str ');'];

  eval(iricall_str)
   
  height_axis = start_height+(0:(num_heights-1))*height_step;

  % Set limits for interpolation
  bottom_transition = 120;        % always start from 120 km, 
				  % ignore top 20 km of FIRI model
  if (iono_extra(4) > 160)
      % if F1 present, interpolate to 5 km below F1 peak
      top_transition = iono_extra(4) - 5;
  else
      % if no F1, interpolate to halfway between 120 km and F2 peak
      top_transition = (bottom_transition + iono_extra(2))/2;
  end

  if (bottom_transition > top_transition)
    bottom_transition = top_transition - 10; 
  end

  firi_transition = ...
      (height_axis > bottom_transition) & (height_axis < top_transition);
  top_idx = find(height_axis > top_transition, 5, 'first');
  bottom_idx = find(height_axis <= bottom_transition, 5, 'last');

  iono(1,firi_transition) = exp(interp1(height_axis([top_idx bottom_idx]), ...
     log(iono(1,[top_idx bottom_idx])), height_axis(firi_transition), 'pchip'));

  valid = iono(1,:) > 0;
  bottom = find(valid, 5, 'first');
  iono(1,~valid) = exp(interp1([40 height_axis(bottom)], ...
      log([1e-10 iono(1,bottom)]), height_axis(~valid), 'linear', 'extrap'));
  iono(1,height_axis < 40) = 0;
end
