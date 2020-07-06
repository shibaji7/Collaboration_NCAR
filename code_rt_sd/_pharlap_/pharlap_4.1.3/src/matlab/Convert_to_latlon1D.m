function [ elats, elons ] = Convert_to_latlon1D(lat, lon, b, D)
% b = bearing in degree
% lat, lon = starting lat & lon in degree
% D = distance in km

dx = D*sin(deg2rad(b));
dy = D*cos(deg2rad(b));

d_lats = dx*1000/(111320*cos(deg2rad(lat)));
d_lons = dy*1000/110540;


elats = lat + abs(d_lats);
elons = lon - abs(d_lons);
end