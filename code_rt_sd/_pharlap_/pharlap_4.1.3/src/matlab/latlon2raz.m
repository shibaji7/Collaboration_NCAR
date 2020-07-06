% Name: 
%   latlon2raz.m
% 
% Purpose:
%   Converts latitude and longitude to polar coordinates in range/azimuth.
%
% Calling sequence:
%   1.  [range, azimuth] = latlon2raz(point_lat, point_lng, origin_lat, ...
%  	      origin_lng)
%   2.  [range, azimuth] = latlon2raz(point_lat, point_lng, origin_lat, ...
%  	      origin_lng, geoid)
%   3.  [range, azimuth] = latlon2raz(point_lat, point_lng, origin_lat, ...
%  	      origin_lng, 'wgs84')
%   4.  [range, azimuth] = latlon2raz(point_lat, point_lng, origin_lat, ...
%  	      origin_lng, 'spherical')
%   5.  [range, azimuth] = latlon2raz(point_lat, point_lng, origin_lat, ...
%  	      origin_lng, 'spherical', re)
%     
% Inputs:
%   point_lat  - latitude of point (array)
%   point_lng  - longitude of point (array)
%    
%   origin_lat - the latitude and longitude of the origin for
%   origin_lng      which the range, azimuth data should refer (scalar).
%  
% Optional Inputs
%   geoid      - The type of ellipsoid to use. If not specfied then a 
%                spherical Earth is assumed. The following values are valid:
%
%     'spherical'   - Spherical Earth with radius = 6378137 m (default) or 
%                     user specified input (see calling sequence example 5.)
%     'wgs84'       - WGS84 
%     'wgs1972'     - WGS 1972
%     'grs80'       - Geodetic Reference System 1980 (GRS 80)
%     'airy'        - Airy
%     'aust_nat'    - Australian National
%     'bessel_1841' - Bessel 1841 
%     'clarke_1866' - Clarke 1866
%     'clarke_1880' - Clarke 1880
%
% Outputs:
%   range   - range as distance from (origin_lat, origin_lng) in meters
%   azimuth - azimuth as true bearings (degrees east of north) from origin
%	
% Notes: 
%   1. The range will never exceed half the circumference of the earth.
%      instead the returned azimuth will be shifted by 180 degrees.
%   2. Uses T. Vincenty's iterative method for the non-spherical geoids. This 
%      should be good over distances from a few cm to nearly 20,000 km with
%      millimeter accuracy.
%   3. For non-spherical geoids the solution will fail to converge near the
%      antipodal point (within 0.01 deg in latitude or 0.1 deg in
%      longitude). For this case the precision will be reduced.
% 
% Modification History:
% 20/07/2006 M. A. Cervera 
%   Initial version - spherical earth only
%
% 26/06/2007 M. A. Cervera
%   Added various ellipsoids: WGS84, WGS 1972, GRS 80, Airy, Australian
%   National, Bessel 1841, Clarke 1866, Clarke 1880. Uses T. Vincenty's
%   method. Converted to matlab from FORTRAN code at 
%   ftp://www.ngs.noaa.gov/pub/pcsoft/for_inv.3d/source/forward.for
%
% 25/10/2012 M. A. Cervera
%   User can now specify the radius of the Earth for the spherical case
%
% 22/03/2013 M. A. Cervera
%   Now returns a warning if near the antipodal point for non-spherical geoids.
%

function  [range, azimuth] = latlon2raz(point_lat, point_lng, origin_lat, ...
      origin_lng, vargin, vargin2)
  
%
% general constants
%
dtor = pi/180;
radeg= 180/pi;


%
% ellipsoid constants 
%
if ~exist('vargin', 'var') vargin = 'spherical'; end

switch lower(vargin)
  case {'wgs84'}
    % WGS84 (www.wgs84.com)
    a = 6378137.0;                  % semi-major axis (m)
    f = 1.0 / 298.257223563;	    % Flattening factor.
 
  case {'airy'}
    % Airy
    a = 6377563.396;
    f = 1.0 / 299.3249646;     
    
  case {'aust_nat'}
    % Australian National
    a = 6378160.0;    
    f = 1.0 / 298.25;          
  
  case {'bessel_1841'}
    % Bessel 1841
    a = 6377397.155;   
    f = 1.0 / 299.1528128;     

  case {'clarke_1866'}
    % Clarke 1866
    a = 6378206.400;   
    f = 1.0 / 294.9786982; 

  case {'clarke_1880'}
    % Clarke 1880
    a = 6378249.145;	
    f = 1.0 / 293.465000;	
  
  case {'grs80'}
    % Geodetic Reference System 1980 (GRS 80)
    a = 6378137.0;      
    f = 1.0 / 298.257222101;   

  case {'wgs1972'}
    % WGS 1972
    a = 6378135.0;      
    f = 1.0 / 298.26;    

  case {'spherical', ''}
    % spherical Earth
    if ~exist('vargin2')
      re = 6378137;        % mean radius of Earth in m
    else
      re = vargin2;
    end
    
  otherwise
     error('unknown geoid')
     return
     
end


%
% do the calculations
%
switch lower(vargin)
  case {'spherical', ''}       % spherical Earth
    bore_lat = origin_lat;    
    bore_long = 180.0;    

    lat = point_lat;
    long = point_lng + 180 - origin_lng;   % as we want longtitude line between
					   % origin and bore to be zero
    chi      = dtor * (90.0 - bore_lat);
    rlat     = dtor * lat;
    rlon     = dtor * long;
    coschi   = cos(chi);
    sinchi   = sin(chi);
    coslat   = cos(rlat);
    coslon   = cos(rlon);
    sinlat   = sin(rlat);
    sinlon   = sin(rlon);

    y        = sinlon.*coslat;
    x        = coslon.*coslat.*coschi + sinchi.*sinlat;

    geo_lat  = asin(coschi.*sinlat - coslon.*coslat.*sinchi) / dtor;
    geo_long = atan2(y, x) / dtor + bore_long;

    % due south is 0 degrees long., +ve anticlockwise,
    % due north is 0 degrees long., +ve clockwise, = true bearings
    azimuth = 180.0 - geo_long;           
    azimuth = mod((azimuth + 360), 360);   % force azimuth to 0 - 360 degrees

    range = abs((90.0 - geo_lat) * dtor * re); 

  otherwise    % ellipsoid
    %---------  Initial values  --------------------
    epsi = 0.5e-13;			% Tolerence.
    glat1 = origin_lat * dtor;
    glon1 = origin_lng * dtor;
    glat2 = point_lat * dtor;
    glon2 = point_lng * dtor;

    r = 1 - f;
    tu1 = r .* sin(glat1) ./ cos(glat1);
    tu2 = r .* sin(glat2) ./ cos(glat2);
    cu1 = 1 ./ sqrt(1 + tu1.^2);
    su1 = cu1 .* tu1;
    cu2 = 1 ./ sqrt(1 + tu2.^2);
    s = cu1 .* cu2;
    baz = s .* tu2;
    faz = baz .* tu1;
    x = glon2 - glon1;

    %-----  Iterate  ---------------
    d = 1e10;
    count = 0;
    while (max(abs(d-x)) > epsi & count < 20)
      count = count+1;
      sx = sin(x);
      cx = cos(x);
      tu1 = cu2 .* sx;
      tu2 = baz - su1.*cu2.*cx;
      sy = sqrt(tu1.^2 + tu2.^2);
      cy = s.*cx + faz;
      y = atan2(sy, cy);
      sa = s.*sx ./ sy;
      c2a = 1 - sa.^2;
      cz = 2 * faz;
      w = find(c2a > 0);
      cnt = length(w);
      if cnt > 0 cz(w) = cy(w) - cz(w) ./ c2a(w); end
      e = 2 * cz.^2 - 1;
      c = ((-3 * c2a + 4) .* f + 4) .* c2a .* f/16;
      d = x;
      x = ((e.*cy.*c + cz) .* sy .* c + y) .* sa;
      x = (1 - c).*x.*f + glon2 - glon1;
    end  
    if count == 20
      warning(['Near the antipodal point - solution failed to converge and ' ...
	       'has reduced precision'])
    end
    
    %------  Finish up  ----------------
    faz = atan2(tu1, tu2);
    baz = atan2(cu1.*sx, baz.*cx - su1.*cu2) + pi;
    x = sqrt((1 ./ r.^2 - 1).*c2a + 1) + 1;
    x = (x - 2) ./ x;
    c = 1 - x;
    c = (x.^2 ./ 4 + 1) ./ c;
    d = (0.375.*x.^2 - 1).*x;
    x = e .* cy;
    s = 1 - e.^2;
    range = ((((4 * sy.^2 - 3) .* s.*cz.*d/6 - x) .* d/4 + cz) .* sy.*d + y) .* ...
	c.*a.*r;

    azimuth = faz ./ dtor;
    w = find(azimuth < 0);
    cnt = length(w);
    if cnt > 0 azimuth(w) = azimuth(w) + 360; end
 
end
