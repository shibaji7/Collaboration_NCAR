%
% Name:
%   noise_ccir.m
%
% Purpose / description:
%   Calculate median atmospheric noise energy at a particular location for a set
%   of UTs and freqs. Based on ITU-R Recommendation ITU-R P.372-11 for radio
%   noise. Galactic and man-made noise are not included - they must be
%   calcualated separately (see Notes below).
%
% Calling Sequence: 
%   atmos_noise = noise_ccir(UT_hour, freqs_mhz, lat, lon, month)
%
% Inputs:
%   UT_hour    -  universal time in hours (vector of length num_times)
%   freqs_mhz  -  model frequency in MHz (vector of length num_freqs)
%   lat        -  latitude in deg (scalar)
%   lon        -  longitude in deg (scalar)
%   month      -  month of year number (scalar)
%
% Outputs - all arrays of size (num_freqs, num_times):
%   atmos_noise -  Fam - 204, Median Level atmospheric noise
%   dla         -  Ratio of lower decile to median value Fam
%   dua         -  Ratio of upper decile to median valuw Fam
%   sla         -  Standard deviation of values of dla
%   sma         -  Standard deviation of values of Fam
%   sua         -  Standard deviation of values of dua
%   vd          -  Voltage deviation
%   svd         -  Standard deviation of Voltage deviation 
%
% References:
%   International Telecommunication Union Recommendation ITU-R P.372-10 
%   (10/2009), Radio Noise. A copy of this may be found in the PHaRLAP 
%   dat/ccir_noise directory. 
%
% Notes:
%   To calculate the total noise at a site the galactic and man-made noise
%   must be added. These may be specified by ITU-R SG3 NOISBW Recommendation
%   ITU-R P.372-11 (see page 12) or by an alternative source. The ITU
%   documentation may be found in the PDF file :
%        R-REC-P.372-11-201309-I!!PDF-E.pdf 
%   located in the folder :
%        <pharlap>/dat/ccir_noise/
%   The following code shows an example where the ITU recommendation has be 
%   use to calculate the total noise.
%
% >> noise_atmos = noise_ccir(UT_hour, radar_freq, radar_lat, radar_lon, month);
% >> noise_galactic = 52 - 23*log10(radar_freq) - 204;     % ITU galactic noise
% >> noise_man = 53.6 - 28.6*log10(radar_freq) - 204;      % ITU quiet rural 
% >> noise_lin = 10.^(noise_atmos/10) + 10.^(noise_galactic/10) + ...
% >>  	         10.^(noise_man/10);
% >> noise_dBW_perHz = 10*log10(noise_lin);
%
%
% Modification History:
%   June 2011, Mike Turley, Defence Science & Technology Organisation - Author
%
%   28/02/2012 M.A.Cervera  
%      Minor modifications apropos integration into PHaRLAP
%
%   03/07/2015 M.A.Cervera  
%      Updated noise recommendation documentation to ITU-R P.372-11
%

function [atmos_noise, dla, dua, sla, sma, sua, vd, svd] = noise_ccir (uts, freqs_mhz, lat, lon, month)
  
  % Constants
  month2season_south = [3,3,4,4,4,1,1,1,2,2,2,3];
  month2season_north = [1,1,2,2,2,3,3,3,4,4,4,1];
  location = 'unknown';
 
  % Size of output arrays
  num_times = length(uts);
  num_freqs = length(freqs_mhz);
  
  % Output array declarations
  atmos_noise = zeros(num_freqs, num_times);  
  dla = zeros(num_freqs, num_times);   
  dua = zeros(num_freqs, num_times);    
  sla = zeros(num_freqs, num_times);    
  sma = zeros(num_freqs, num_times);    
  sua = zeros(num_freqs, num_times);   
  vd = zeros(num_freqs, num_times);     
  svd = zeros(num_freqs, num_times);   
  
  % Make sure this lon belongs to this day (0-360)
  lon = mod(lon - 1, 360) + 1;
      
  % Loop over time of day
  for itime = 1:num_times
    
    local_time = mod(uts(itime) + lon / 360 * 24 - 1, 24) + 1; 
    if lat < 0
      season = month2season_south(month);
    else
      season = month2season_north(month);
    end
    
    for ifreq = 1:num_freqs
      
      % Call matlab implementation of the NOISEDAT main program NOISBW
      apply_correction = false;
      noisbw = radio_noise(season, lat, lon, location, freqs_mhz(ifreq), ...
	                   local_time, false, apply_correction);
      
      % Save into arrays
      atmos_noise(ifreq,itime) = [noisbw.atmos];
      dla(ifreq,itime) = [noisbw.dlas];
      dua(ifreq,itime) = [noisbw.duas];
      sla(ifreq,itime) = [noisbw.slas];
      sma(ifreq,itime) = [noisbw.smas];
      sua(ifreq,itime) = [noisbw.suas];
      vd(ifreq,itime)  = [noisbw.vds];
      svd(ifreq,itime) = [noisbw.svds];
      
    end % ifreq
  end % itime
  
  % Done

return