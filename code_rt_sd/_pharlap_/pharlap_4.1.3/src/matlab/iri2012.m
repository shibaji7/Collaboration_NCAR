%
% Name :
%   iri2012.m
%
% Purpose :
%   Matlab wrapper to the IRI-2012 fortan based empirical model ionosphere. 
%   For (1) a given user specified R12 index (i.e. ionospheric condition), or 
%   (2) a given epoch, or (3) user specified ionospheric layers, it returns
%   modeled electron density, electron temperature, ion temperature, and ion
%   composition (O+, H+, He+, NO+, O+2) height profiles, ion drift and F1
%   probability along with various other parameters including foE, hmE, foF1,
%   hmF1, foF2, hmF2. Storms may be switched on or off.
%   NB: SEE NOTES SECTION BELOW FOR IMPORTANT INFORMATION.  
%
% Calling sequence :
%   1. [iono, iono_extra] = iri2012(lat, lon, R12, UT);
%
%       Ionospheric profiles are not generated, ionospheric parameters (foF2,
%       hmf2 etc.) only are returned.
%
%   2. [iono, iono_extra] = iri2012(lat, lon, R12, UT, ht_start, ht_step ...
%                                    num_hts);
%       Ionospheric profiles generated and returned in output iono along with
%       the iononospheric parameters.
%
%   3. [iono, iono_extra] = iri2012(lat, lon, R12, UT, ht_start, ht_step ...
%                                    num_hts, B0B1_model);
%       Ionospheric profiles generated and returned in output iono along with
%       the iononospheric parameters. User selected model for the ionospheric 
%       profile.
%
%   4. [iono, iono_extra] = iri2012(lat, lon, R12, UT, ht_start, ht_step ...
%                                    num_hts, B0B1_model, D_model);
%       Ionospheric profiles generated and returned in output iono along with
%       the iononospheric parameters. User selected model for the ionospheric 
%       profile. User selected model for D Layer.
%
%   5. [iono, iono_extra] = iri2012(lat, lon, R12, UT, ht_start, ht_step ...
%                                    num_hts, B0B1_model, D_model, ...
%                                    iono_layer_parms);
%       Ionospheric profiles generated and returned in output iono along with
%       the iononospheric parameters. User selected model for the ionospheric
%       profile. User selected model for D Layer. The E, F1 & F2 layer
%       parameters are user specified in the input iono_layer_parms. 
%
% Inputs :
%   lat      - geographic latitude of point (degrees)  
%   lon      - geographic longitude of point (degrees)
%   R12      - scalar or 2 element array :
%    scalar:
%      R12 = 1 - 200 :  IRI2012 is called with R12 input as the user specified
%                       yearly smoothed monthly median sunspot number. The
%                       storm model is turned off
%      R12 = -1      :  IRI2012 is called with ionospheric conditions read
%                       from file (ig_rz.dat) based on input epoch (UT) and
%                       may be historical or projected conditions (dependent
%                       on the epoch). The storm model is turned off. See
%                       also notes below.
%     array: [R12 storm_flag]
%       1st element is R12 index as above
%       2nd element is a flag to turn storm model on (storm_flag = 1) off (= 0)
%
%   UT       - 5x1 array containing UTC date and time - year, month, day, 
%              hour, minute. The year must be in the range 1958 - 2018
%
% Optional Inputs:
%   ht_start - start height for ionospheric profile generation (km)
%   ht_step  - height step size for ionospheric profile generation (km)
%   num_hts  - number of heights for ionospheric profile generation
%              must be >= 2 and <= 1000
%
%   B0B1_model - specifies the model to use for the bottomside ionospheric
%                profile. The model defines how IRI2012 determines the B0
%                (thickness) and B1 (shape) parameters. If B0B1_model is not 
%                specified then the default is Bil-2000 (B0B1_model = 2)
%       B0B1_model = 1 : ABT-2009  (Adv. Space Res., Vol 43, 1825-1834, 2009)
%       B0B1_model = 2 : Bil-2000  (Adv. Space Res., Vol 25, 89-96, 2000)
%       B0B1_model = 3 : Gul-1987  (Adv. Space Res., Vol 7, 39-48, 1987)
%
%   D_model - specifies the model to use for the D Layer. If D_model is not 
%             specified then the default is IRI-1990 (D_model = 1)
%       D_model = 1 : IRI-1990
%       D_model = 2 : FT-2001 and DRS-1995 - Friedrich and Torkar's FIRI model
%                     for the lower ionosphere (Friedrich and Torkar, 
%  	              J. Geophys Res. 106 (A10), 21409Ð21418, 2001)
%                     Danilov's et al. D-region model (Danilov et al, Adv. 
%		      Space Res., vol 15, 165, 1995) is also returned in 
%                     the output iono(14, :)
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
%   iono  -  Array of 14 output quantities x num_hts heights. NB, ionosperic
%            profiles are not computed below 65 km or above 2000 km. Values
%            of -1 are returned if these heights are requested and also for
%            those heights where a valid value is unable to be calculated. If
%            the optional inputs  ht_start and ht_step are not supplied then
%            the profiles are not calculated.
%       iono(1, :) = electron number density (m^-3). If D_model = 2 then
%                    Ne is given by the Friedrich and Torkar (2001) FIRI model 
%                    in the height interval 60 - 140 km
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
%       iono(12, :) = NOT USED (value = -1)
%       iono(13, :) = NOT USED (value = -1)
%       iono(14,  :  )   if D_model = 1 then not used (value = -1), otherwise
%                1:11) = standard IRI Ne for 60, 65, ... 110km 
%               12:22) = Friedrich and Torkar FIRI model at these heights 
%               23:33) = standard Danilov (SW=0, WA=0) at these heights 
%               34:44) = Danilov for minor Stratospheric Warming (SW=0.5) 
%               45:55) = Danilov for major Stratospheric Warming (SW=1) 
%               56:66) = Danilov weak Winter Anomaly (WA=0.5) conditions
%               67:77) = Danilov strong Winter Anomaly (WA=1) conditions
%               78:  ) = NOT USED (value = -1)
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
%   2.3 IRI indices files ig_rz.dat and apf107.dat updated 04-May-2016
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
%   V1.0  M. A. Cervera  05/09/2012
%      Initial version
%   
%   V1.1  M. A. Cervera  16/12/2013
%      Now accepts user input of foF2, hmF2, foF1, hmF1, foE and hmE
%      if required
%
%   V1.2  M. A. Cervera  25/06/2015
%      Updated to latest version of IRI2012. Users can now specify if the
%      storm model is to be turned on or off. The model for the ionospheric 
%      profile may now be selected. Te model for the D Layer may now be
%      seleceted.
%

% This a Matlab help file only. The actual programme is a mex wrapper
% (iri2012_matlab_wrapper.for) to the Fortran code (irisub.for).
