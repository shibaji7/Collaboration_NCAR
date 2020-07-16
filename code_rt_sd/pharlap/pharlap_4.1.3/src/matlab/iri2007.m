%
% Name :
%   iri2007.m
%
% Purpose :
%   Matlab wrapper to the IRI-2007 fortan based empirical model ionosphere. 
%   For (1) a given user specified ionospheric condition, or (2) a given epoch
%   it returns modeled electron density, electron temperature, ion temperature,
%   and ion composition (O+, H+, He+, NO+, O+2) height profiles, ion drift
%   and F1 probability along with various other parameters including foE, hmE, 
%   foF1, hmF1, foF2, hmF2. Storms are not modelled for case (1) and may be 
%   switched on or off for case (2). NB: SEE NOTES BELOW FOR IMPORTANT 
%   INFORMATION.
%
% Calling sequence :
%    1. [iono, iono_extra] = iri2007(lat, lon, R12, UT, ht_start, ht_step);
%    2. [iono, iono_extra] = iri2007(lat, lon, R12, UT);
%
% Inputs :
%   lat      - geographic latitude of point (degrees)  
%   lon      - geographic longitude of point (degrees)
%   R12      - depends on value as follows:
%      R12 = 1 - 200 :  IRI2007 is called with R12 input as the user specified
%                       yearly smoothed monthly median sunspot number. The
%                       storm model is turned off
%      R12 = -1      :  IRI2007 is called with ionospheric conditions read
%                       from file (ig_rz.dat) based on input epoch (UT) and
%                       may be historical or projected conditions (dependent
%                       on the epoch). See notes below.
%      R12 = -2      :  as for R12 = -1 but storm model is turned off 
%
%   UT       - 5x1 array containing UTC date and time - year, month, day, 
%                hour, minute
%
% Optional Inputs:
%   NOTE: If the optional inputs are to be specified then they all must supplied
%   ht_start - start height for ionospheric profile generation (km)
%   ht_step  - height step size for ionospheric profile generation (km)
%
% Outputs :
%   iono  -  Array of 11 output parameters x 100 heights. NB, ionosperic
%            profiles are not computed below 60 km. Values of -1 are returned
%            if these heights are requested and also for those heights where
%            a valid value is unable to be calculated. If the optional inputs
%            ht_start and ht_step are not supplied then the profiles are not
%            calculated and -1 is returned for all elements.
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
%      iono_extra(43) = daynr
%      iono_extra(44) = equatorial vertical ion drift in m/s
%      iono_extra(45) = foF2_storm/foF2_quiet
%      iono_extra(46) = F1 probability without L condition 
%      iono_extra(47) = F1 probability with L condition incl.
%      iono_extra(48 - 50) = Unused
%
%  Notes :
%   1. Notes for IRI2007 called using specified input ionospheric conditions:
%   1.1 If the ionospheric conditions are controlled by the matlab input R12
%      index then the input year (in the UT array) has no effect on the solar
%      conditions. However, it does have a small affect on the magnetic field
%      model which could be important.
%
%   1.2 User defined IG12 and F10.7 are derived from R12 using the following 
%      empirical formulas :
%            F107 = 63.7 + 0.728*R12 + 0.00089*R12^2
%            IG12 = -12.349154 + R12 * (1.4683266 - R12 * 2.67690893e-03);
%      These derived values for IG12 and F10.7 are input into IRI-2007
%
%   2. Notes for IRI2007 called using specified input epoch:
%   2.1 IRI2007 uses solar indices tabled in ig_rz.dat (which is supplied with
%      IRI2007). This file contains R12 and IG12 from Jan 1958 to Dec 2011. If 
%      the input UT is outside this range then an error is returned. R12 is
%      used to model the height of the F2 layer and IG12 its strength.
%
%   2.2 The computation of the yearly-running mean for month M requires the
%      indices for the six months preceeding M and the six months following 
%      M (month: M-6, ..., M+6). To calculate the current running mean one 
%      therefore requires predictions of the index for the next six months. 
%      Starting from six months before the UPDATE DATE (listed at the top of 
%      the file ig_rz.dat and below at point 3.) and onward the indices are 
%      therefore based on indices predictions.
%
%   2.3 ig_rz.dat updated 25-Aug-2008
%
%   2.4 The solar activity parameter F10.7 (daily) and magnetic activity Ap index
%      (3-hourly) used by IRI2007 are tabled in ap.dat from 1 Jan 1960 to 30 Jun 
%      2008. If UT is outside this range then the storm model (which relies on
%      Ap) is switched off and a monthly median F10.7 (calculated from R12 using
%      the empirical formula F107 = 63.7 + 0.728*R12 + 0.00089*R12^2) is used in
%      place of the daily F10.7. F10.7 is used for the topside electron 
%      temprature model (Intercosmos) and equatorial vertical ion drift model.
%
%   3. This mex file drives IRI-2007 with the following options
%      5.1  B0   - Table option
%      5.2  foF2 - URSI
%      5.3  Ni   - DS-95 and TTS-03
%      5.4  Ne   - f10.7 unlimited
%      5.5  IGRF magnetic field model
%      5.6  F1 probability model
%      5.7  standard F1
%      5.8  ion drift computed
%      5.9  Te_topside - Intercosmos
%      5.10 D-Region model - IRI-95
%      5.11 spread-F probability - not computed
%      5.12 NeQuick topside model
%      5.13 R12, IG12 and F107D are user defined for R12 = 1-200, otherwise 
%           read from file
%      5.14 Storm model is turned on for R12 = -1, otherwise off
%       
%
% Further information :
%   http://iri.gsfc.nasa.gov/
%
% Modification History:
%   V1.0  M.A. Cervera  28/10/2008
%
%   V1.1  M.A. Cervera  09/12/2009
%      ht_start and ht_step are now optional inputs, if they are not supplied
%      then the profiles are not calculated which may speed up the call to IRI
%

% This a Matlab help file only. The actual programme is a mex wrapper
% (iri2007_matlab_wrapper.for) to the Fortran code (irisub.for).
