%
% Name :
%   iri2016.m
%
% Purpose :
%   Matlab wrapper to the IRI-2016 fortan based empirical model ionosphere. 
%   For (1) a given user specified R12 index (i.e. ionospheric condition), or 
%   (2) a given epoch, or (3) user specified ionospheric layers, it returns
%   modeled electron density, electron temperature, ion temperature, and ion
%   composition (O+, H+, He+, NO+, O+2) height profiles, ion drift and F1
%   probability along with various other parameters including foE, hmE, foF1,
%   hmF1, foF2, hmF2. Storms may be switched on or off.
%   NB: SEE NOTES SECTION BELOW FOR IMPORTANT INFORMATION.  
%
% Calling sequence :
%   1. [iono, iono_extra] = iri2016(lat, lon, R12, UT);
%
%       Ionospheric profiles are not generated, ionospheric parameters (foF2,
%       hmf2 etc.) only are returned.
%
%   2. [iono, iono_extra] = iri2016(lat, lon, R12, UT, ht_start, ht_step ...
%                                    num_hts);
%       Ionospheric profiles generated and returned in output iono along with
%       the iononospheric parameters.
%
%   3. [iono, iono_extra] = iri2016(lat, lon, R12, UT, ht_start, ht_step ...
%                                    num_hts, iri_options)
%       Ionospheric profiles generated and returned in output iono along with
%       the iononospheric parameters. Defauls IRI behaviour overridden by 
%       iri_options.  
%
% Inputs :
%   lat      - geographic latitude of point (degrees)  
%   lon      - geographic longitude of point (degrees)
%   R12      - scalar R12 index
%      R12 = 1 - 200 :  IRI2016 is called with R12 (Zurich V1.0) input as the
%                       user specified yearly smoothed monthly median sunspot
%                       number. The foF2 storm model will be turned off
%                       regardless of the setting of the optional input
%                       iri_options.foF2_storm (see below). 
%      R12 = -1      :  IRI2016 is called with ionospheric conditions (R12,
%                       IG12, and F10.7) read from file (ig_rz.dat) based on
%                       input epoch (UT) and may be historical or projected
%                       conditions (dependent on the epoch). The foF2 storm
%                       model defaults to "on" for this case (but may be
%                       overridden). See the optional input
%                       iri_options.foF2_storm below.   
%
%      NB: IRI also requires the IG12 and F10.7 ionospheric parameters to be
%          supplied. For R12 = -1 these will be read from file (as for R12). 
%          For a user supplied R12 they are calculated in the mex wrapper
%          from R12 and input into IRI. The formulae used are:
%                 F107 = 63.75 + R12 * (0.728 + R12*0.00089)
%                     (see Davies, "Ionospheric Radio", 1990, pp442)  
%                 IG12 = -12.349154 + R12 * (1.4683266 - R12*2.67690893e-03)
%                     (see irisub.for, line 913)
%
%   UT       - 5x1 array containing UTC date and time - year, month, day, 
%              hour, minute. The year must be in the range 1958 - 2018
%
% Optional Inputs:
%   ht_start          - start height for ionospheric profile generation (km)
%   ht_step           - height step size for ionospheric profile generation (km)
%   num_hts           - number of heights for ionospheric profile generation
%                       must be >= 2 and <= 1000
%
%   iri_options - structure containing options to control IRI. If iri_options
%                 is not used or if a field has not been specfied or has an
%                 invalid value then the default value is  used. Valid fields
%                 are : 
%       .iri_messages - turns IRI messages on or off
%                         'off' - messages off (default)
%                         'on'  - messages on
%
%       .foF2 - user input foF2. Leave undefined to use IRI model.
%               Note 1. range of values : 0.1Mz < foE < foF1 < foF2 < 100MHz
%                    2. If defined then .foF2_storm will be ignored 
%
%       .hmF2 - user input hmF2. Leave undefined to use IRI model.
%               Note range of values : 50km < hmE < hmF1 < hmF2 < 1000km
%
%       .foF1 - user input foF2. Leave undefined to use IRI model.
%               Note range of values : 0.1Mz < foE < foF1 < foF2 < 100MHz
%
%       .hmF1 - user input hmF1. Leave undefined to use IRI model.
%               Note range of values : 50km < hmE < hmF1 < hmF2 < 1000km
%
%       .foE  - user input foE. Leave undefined to use IRI model.
%               Note range of values : 0.1Mz < foE < foF1 < foF2 < 100MHz
%
%       .hmE  - user input hmE. Leave undefined to use IRI model.
%               Note range of values : 50km < hmE < hmF1 < hmF2 < 1000km
%
%       .foF2_coeffs - specifies which coefficients to use for foF2 model.
%                      Valid values are:
%                        'URSI' - (default)
%                        'CCIR'
%
%       .Ni_model - specifies models to use for ion density profile. Valid
%                   values are 
%                     'RBV-2010 & TTS-2005' (default)
%                     'DS-1995 & DY-1985'
%
%       .Te_profile - specifies whether to use standard Te or Te/Ne correlation
%                  'Te/Ne correlation' - Te calculated using Te/Ne correlation 
%                  'standard'          - standard Te calculation (default)
%       
%       .Te_topside_model - specifies Te topside model
%                       'TBT-2012' - (default)
%                       'Bil-1985'
%
%       .Te_PF107_dependance - specifies whether to use PF10.7 dependance for 
%                              the Te model
%                     'off' - no PF10.7 dependance
%                     'on'  - with PF10.7 dependance (default)
%
%       .Ne_tops_limited - specifies whether or not f10.7 is limited for the
%                          purpose of calculating the topside electron density
%                            'f10.7 unlimited' - f10.7 unlimited
%                            'f10.7 limited'   - f10.7 limited to 188 (default)
%
%       .Ne_profile_calc - specifies whether to use standard Ne calculation or
%                          Ne using Lay-function formalism.
%                   'Lay-function' - Ne calculated using Lay-function formalism
%                   'standard'     - standard Ne calculation (default)
%
%       .Ne_B0B1_model - string which specifies the model to use for the
%                     bottomside ionospheric profile. The model defines how
%                     IRI2016 determines the B0 (thickness) and B1 (shape)
%                     parameters. The default is ABT-2009. Valid values are:
%                       'ABT-2009'  (Adv. Space Res., Vol 43, 1825-1834, 2009) 
%                       'Bil-2000'  (Adv. Space Res., Vol 25, 89-96, 2000)
%                       'Gul-1987'  (Adv. Space Res., Vol 7, 39-48, 1987)
%
%       .Ne_topside_model - string which specifies the model to use for the
%                           topside ionospheric profile. Valid values are:
%                      'IRI-2001' - topside model from IRI2001
%                      'IRI-2001 corrected' - corrected IRI2001 topside model
%                      'NeQuick' - the NeQuick topside model (default)
%  
%       .F1_model - specifies model to use for F1 layer. See Scotto et al., 
%                   Adv. Space Res., Vol 20, Number 9, 1773-1775, 1997.
%                     'Scotto-1997 no L' - Scotto without L condition (default)
%                     'Scotto-1997 with L' - Scotto with L condition
%                     'solar zenith' - critical solar zenith angle (old IRI95)
%
%       .D_model - specifies the model to use for the D Layer. The default is
%                  IRI-1990. Valid values are: 
%                    'IRI-1990' - Model from IRI90 (default)
%                    'FT-2001'  - Friedrich and Torkar's FIRI model for the 
%                                 lower ionosphere (Friedrich and Torkar,
%                                 J. Geophys Res. 106 (A10), 21409Ã?21418,
%                                 2001). Danilov's et al. D-region model
%                                 (Danilov et al,  Adv. Space Res., vol 15,
%                                 165, 1995) is also returned in the output
%                                 iono(14, :) 
%
%       .hmF2_model - specifies the model to use for hmf2. Default is AMTB.
%                     Valid values are: 
%                       'AMTB' (default)
%                       'Shubin-COSMIC'
%                       'M3000F2'
%
%       .foF2_storm - specifies whether or not to have the foF2 storm model on.
%                     NB: If either R12 has been supplied or foF2 has been user 
%                     input (i.e. field .foF2 is set) then this will be 
%                     ignored as it is no longer relevant and the foF2 storm
%                     model is turned off. 
%                       'off' - no storm model 
%                       'on'  - storm model on (default)
%
%       .hmF2_storm - specifies whether or not to have foF2 storm model on
%                     for hmF2 model
%                       'off' - no storm model (default)
%                       'on'  - storm model on
%
%       .foE_storm - specifies whether or not to have foE storm model on
%                       'off' - no storm model (default)
%                       'on'  - storm model on
%
%       .topside_storm - specifies whether or not to have foF2 storm model on
%                        for the topside model
%                          'off' - no storm model (default)
%                          'on'  - storm model on
%       
%       .auroral_boundary_model - specifies whether to have auroral boundary
%                                 model on or off
%                          'off' - auroral boundary model is off (default)
%                          'on'  - auroral boundary model is on
%       
%       .covington - method for calculating Covington Index. Valid values are:
%                      'F10.7_12' - (default)
%                      'IG12' - used by IRI before Oct 2015
%
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
%       iono(14,  :  )   if iono_options.D_model = 'IRI-1990' then not used 
%                        (value = -1), otherwise :
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
%      ionohh_extra(53) = invdip/degree      iono_extra(54) = MLT
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
%   1. Notes for IRI2016 called using specified input ionospheric conditions:
%   1.1 If the ionospheric conditions are controlled by the matlab input R12
%       index then the input year (in the UT array) has a very small effect
%       on  the solar conditions. For example 
%          [iono iono_extra] = iri2016(-25, 135, 70, [2000 1 1 3 0])
%       returns NmF2 = 1.0252e+12 electrons/m^3, whereas
%          [iono iono_extra] = iri2016(-25, 135, 70, [2001 1 1 3 0])
%       returns NmF2 = 1.0260e+12
%
%   1.2 User defined IG12, F10.7 and F10.7_81 (3 solar rotation average of
%       F10.7) required by IRI2016 are derived from R12 using the following
%       empirical formulas : 
%            F107    = 63.7 + 0.728*R12 + 0.00089*R12^2
%            F107_81 = F107
%            IG12 = -12.349154 + R12 * (1.4683266 - R12 * 2.67690893e-03);
%       These derived values for IG12, F10.7 and F10.7_81 are input into 
%       IRI-2016
%
%   2. Notes for IRI2016 called using specified input epoch:
%   2.1 IRI2016 uses solar indices tabled in ig_rz.dat (which is supplied with
%       IRI2016). This file contains R12 and IG12 from Jan 1958 to Dec 2018. If 
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
%       index (3-hourly) used by IRI2016 are tabled in apf107.dat from 1 Jan 
%       1958 to 31 Dec 2014. If UT is outside this range then the storm model
%       (which relies on  Ap) is switched off and a monthly median F10.7
%       (calculated from R12 using the empirical formula F107 = 63.7 +
%       0.728*R12 + 0.00089*R12^2) is used in  place of the daily
%       F10.7. 
%
%   3. This mex file drives IRI-2016 with the following default options
%      unless over-ridden and input by user via the specification of the 
%      optional input iri_options.
%      3.1  Ne computed
%      3.2  Te, Ti computed
%      3.3  Ne & Ni computed
%      3.4  B0,B1 - other models (set by 3.31)
%      3.5  foF2  - URSI 
%      3.6  Ni    - RBV-2010 & TTS-2005
%      3.7  Ne    - f10.7 unlimited
%      3.8  foF2 from model
%      3.9  hmF2 from model
%      3.10 Te    - Standard
%      3.11 Ne    - Standard Profile
%      3.12 Messages to unit 6 (but see 3.34 below)
%      3.13 foF1 from model
%      3.14 hmF1 from model
%      3.15 foE  from model
%      3.16 hmE  from model
%      3.17 Rz12 from file (unless over-ridden and input by user, i.e. R12 > 0)
%      3.18 IGRF magnetic field model
%      3.19 F1 probability model
%      3.20 standard F1
%      3.21 ion drift computed
%      3.22 ion densities in m-3
%      3.23 Te_topside - TBT-2012
%      3.24 D-Region model - IRI-1990
%      3.25 F107D from APF107.DAT (unless over-ridden by user, i.e. R12 > 0)
%      3.26 foF2 with storm model
%      3.27 IG12 from file (unless over-ridden and input by user, i.e. R12 > 0)
%      3.28 spread-F probability - computed
%      3.29 false i.e. topside defined by 3.30 below
%      3.30 NeQuick topside model 
%      3.31 B0,B1 model set to ABT-2009
%      3.32 F10.7_81 from file  (unless over-ridden by user, i.e. R12 > 0)
%      3.33 Auroral boundary model is off
%      3.34 Messages off
%      3.35 no foE storm updating
%      3.36 hmF2 without foF2-storm                
%      3.37 topside without foF2-storm 
%      3.38 turn WRITEs off in IRIFLIP
%      3.39 hmF2 (M3000F2) false - i.e. new models selected
%      3.40 hmF2 AMTB model selected
%      3.41 Use COV=F10.7_12
%      3.42 Te with PF10.7 dependance
%
% Further information :
%   http://iri.gsfc.nasa.gov/
%   http://irimodel.org/
%
% Modification History:
%   V1.0  M. A. Cervera  26/02/2016
%      Initial version
%   
%

% This a Matlab help file only. The actual programme is a mex wrapper
% (iri2016_matlab_wrapper.c) to the Fortran code (irisub.for).
