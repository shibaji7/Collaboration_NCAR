%
% Name:
%   raz2latlon.m
%
% Purpose:
%   Converts the ground range and azimuth (from true North) of point from a 
%   particular bore point or origin on the Earth to latitude and longitude.
%	
% Calling sequence:
%   1.  [geo_lat, geo_long] = raz2latlon(range, azim, bore_lat, bore_long)
%   2.  [geo_lat, geo_long] = raz2latlon(range, azim, bore_lat, bore_long, ...
%                                        geoid)
%   3.  [geo_lat, geo_long] = raz2latlon(point_lat, point_lng, origin_lat, ...
%  	                                 origin_lng, 'wgs84')
%   4.  [geo_lat, geo_long] = raz2latlon(point_lat, point_lng, origin_lat, ...
%  	                                 origin_lng, 'spherical')
%   5.  [geo_lat, geo_long] = raz2latlon(point_lat, point_lng, origin_lat, ...
%  	                                 origin_lng, 'spherical', re)
%
%
% Inputs:	
%   range 	Array	range (m) of point from the bore (or origin). Bore is the
%			"North Pole" in the bore hole system.
%   azim 	Array	azimuth of point from true North with respect to the 
%                       origin (degrees).
%   bore_long	Scalar	longitude of origin (degrees).	
%   bore_lat	Scalar	latitude of origin (degrees).
%
% Optional Inputs
%   geoid       String  The type of ellipsoid to use. If not specfied then a 
%                       spherical Earth is assumed. The following values are
%                       valid:
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
%   geo_lat	Array	latitude of point (degrees).
%   geo_long	Array	longitude of point (degrees).
%
% Notes:
%   1. Uses T. Vincenty's method for the non-spherical geoids. This should be
%      good over distances from a few cm to nearly 20,000 km with millimeter 
%      accuracy.
%
% Modification History:
% 09/02/1996 M. A. Cervera 
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
% 17/12/2013 L. H. Pederick
%   Minor modification for efficiency - added 'var' argument to exist()
%   function calls

function [geo_lat,geo_long] = raz2latlon(range, azim, bore_lat, ...
                                         bore_long, vargin, vargin2)
  
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
    if ~exist('vargin2', 'var')
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
     latbore = pi/2 - range/re;            %lat. of point in bore system
     longbore = dtor*(180 - azim);         %long. of point in bore system

     chi        = dtor*(90.0 - bore_lat);
     coschi     = cos(chi);
     sinchi     = sin(chi);
     coslatb    = cos(latbore);
     coslonb    = cos(longbore);

     sinlatb    = sin(latbore);
     sinlonb    = sin(longbore);

     geo_lat    = radeg*asin(coschi*sinlatb - coslonb.*coslatb*sinchi);
     y          = sinlonb.*coslatb;
     x          = coslonb.*coslatb*coschi + sinchi*sinlatb;

     geo_long   = radeg*atan2(y,x) + bore_long;
  
  
  otherwise    % ellipsoid
    %---------  Initial values  --------------------
    epsi = 0.5e-13;			% Tolerence.
    glon1 = bore_long * dtor;
    glat1 = bore_lat * dtor;
    faz = azim * dtor;

    r = 1 - f;
    tu = repmat(r * sin(glat1) / cos(glat1), size(faz));
    sf = sin(faz);
    cf = cos(faz);
    baz = 0 * faz;		% Want as many baz az faz.
    w = find(cf ~= 0);
    cnt = length(w);
    if (cnt > 0) 
      baz(w) = atan2(tu(w), cf(w)) * 2;
    end

    cu = 1 ./ sqrt(tu.^2 + 1);
    su = tu .* cu;
    sa = cu .* sf;
    c2a = 1 - sa.^2;
    x = sqrt((1 ./ r.^2 - 1) .* c2a + 1) + 1;
    x = (x - 2) ./ x;
    c = 1 - x;
    c = (x.^2 / 4 + 1) ./ c;
    d = (0.375 * x.^2 - 1) .* x;
    tu = range ./ r ./ a ./ c;
    y = tu;

    %------  Iterate  --------------
    while (max(abs(y-c)) > epsi)
      sy = sin(y);
      cy = cos(y);
      cz = cos(baz + y);
      e = 2 * cz.^2 - 1;
      c = y;
      x = e.*cy;
      y = 2*e - 1;
      y = ( (y .* cz .* d .* (4 * sy.^2 - 3)/6 + x) .* d/4 - cz) .* sy .*d + tu;
    end

    %-------  Finish up  -----------
    baz = cu.*cy.*cf - su.*sy;
    c = r * sqrt(sa.^2 + baz.^2);
    d = su.*cy + cu.*sy.*cf;
    glat2 = atan2(d, c);
    c = cu.*cy - su.*sy.*cf;
    x = atan2(sy.*sf, c);
    c = ( (-3*c2a + 4) .* f + 4) .* c2a .* f/16;
    d = ( (e.*cy.*c + cz) .* sy .* c + y) .* sa;
    glon2 = glon1 + x - (1 - c) .* d .* f;
    baz = atan2(sa, baz) + pi;

    geo_long = glon2 ./ dtor;
    geo_lat = glat2 ./ dtor;
    azi2 = baz ./ dtor;
     
 end

