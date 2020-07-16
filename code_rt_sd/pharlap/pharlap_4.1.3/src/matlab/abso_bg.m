%
% Name : 
%   abso_bg.m
%
% Purpose :
%   Calculates absorption loss due to traversing the D-region via George and
%   Bradley (1974).
%
% Calling sequence :
%   absorption = abso_bg(lat, lon, elev, freq, UT, R12, O_mode);
%
% Inputs :
%   lat      - latitude of mid-point of ray (deg)  
%   lon      - longitude of mid-point of ray (deg)
%   elev     - initial elevation of the ray (deg)
%   freq     - radio frequency of the ray (MHz)
%   UT       - universal time, 5x1 array (year, month, day, hour, minute)
%   R12      - R12 index
%   O_mode   - flag to indicate whether to calculate absorption on the O-mode 
%              ray (logical true) or X-mode ray (logical false)
%
% Outputs :
%   absorption - the loss due to D-Region absorption (dB)
%
% Notes :
%   1. The George and Bradley (1974) model of absorption is based around
%      measurements of absorption at vertical incidence and a method of
%      translating these absorption measurements to oblique paths. However,
%      the model is only designed for typical oblique paths and so it is 
%      not expected to be accurate for “non-typical” oblique paths
%      (e.g. Chordal or ducted paths, paths where the ray spends a long
%      time in a cusp, and high-ray modes) or trans-ionospheric paths. For
%      such ray paths the absorption must be calculated along the ray path
%      (see e.g. Pederick and Cervera, 2014).
%
%
% References :
%   1. George, P. L., and P. A. Bradley (1973), Relationship between H.F. 
%      absorption at vertical and oblique incidence, Proc. Inst. Electr. Eng.,
%      120(11), 1355–1361, doi:10.1049/piee.1973.0273.
%
%   2. George, P. L., and P. A. Bradley (1974), A new method of predicting the
%      ionospheric absorption of high frequency waves at oblique incidence,
%      Telecommun. J., 41(5), 307–311. 
%
%   3. Davies, K. (1990), Ionospheric Radio, pages 222-226, The Institution of
%      Engineering and Technology, London, United Kingdom. 
%
%   4. Pederick, L. and M. A. Cervera (2014), Semi-empirical Model for
%      Ionospheric Absorption based on the NRLMSISE-00 absorption model (SiMIAN),
%      Radio Sci., 49, 81--93, doi:10.1002/2013RS005274
%
%

% Modification history:
%   V1.0  M.A. Cervera  11/09/2006
%
%   V1.1  M.A. Cervera  11/09/2009
%     Fixed  error where input lat,lon was at the start of ray rather than the
%     mid point

% This a Matlab help file only. The actual programme is a mex wrapper
% (abso_bg_matlab_wrapper.for) to the Fortran code (abso_bg.for). 
