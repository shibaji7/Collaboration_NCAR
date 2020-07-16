%
% Name:
%   land_sea                                
%
% Purpose:
%   Routine to define land or sea globally. Resolution is 0.1 degrees.    
%                                                                       
% Calling sequence:
%   terrain_type = land_sea(glat, glong)
%
% Inputs: 
%   glat  - geographic latitude  (-90 to 90 degrees)             
%   glong - geographic longitude (0 to 360 degrees measured East from prime
%           meridien)                            
%                                                                       
% Output: 
%   terrain_type - 0 = sea, 1 = land                             
%                                                                       
% Dependencies:
%   None.
%
% Modification history:
%   27/10/2005  V1.0  M. A. Cervera 
%     Initial version.
%
%   06/05/2010  V1.1  D. J. Netherway  
%     Allow arguments to be equal sized arrays
%

function terrain_type = land_sea(glat, glong)

  persistent map_data

  % Obtain the reference data directory fom the relevant environment 
  refdata_dir = getenv('DIR_MODELS_REF_DAT');

  % open the data file and read in the land/sea data if this is the first
  % call to land_sea.m
  if isempty(map_data)
    filename = [refdata_dir '/global_land_Mask_3600_by_1800.dat'];
    fid = fopen(filename, 'r');    

    if (fid == -1) 
      error('%s\n%s\n%s', ...
	  'Missing data file : global_land_Mask_3600_by_1800.dat', ...
	  'Please set the DIR_MODELS_REF_DAT environment variable to point ',...
	  'to the data file directory')
    end

    map_data = fix(fscanf(fid, '%c', [3600, 1800]));
    
    fclose(fid);
  end

  % determine whether the lat/long is land or sea.
  zlong = mod(glong, 360);
  zlat = glat + 90;
  ilat = mod(round(10*zlat), 1800) + 1;
  ilong = mod(round(10*zlong), 3600) + 1;
  
  % terrain_type = map_data(ilong, ilat);

  % Ensure arrays are one dimensional
  ilong = ilong(:);
  ilat = ilat(:);

  nrows = size(map_data, 1);
  
  terrain_type = NaN(size(ilat));
  valid_terain = ~isnan(ilat);
  
  terrain_type(valid_terain) = map_data(ilong(valid_terain) + (ilat(valid_terain)-1)*nrows);
  
  terrain_type = reshape(terrain_type, size(glat));
  
return
end

