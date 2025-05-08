function [xd,yd, d, az] = latlon_to_disaz(lat1,lon1,lat2,lon2)
%  LATLON_TO_DISAZ Haversine formula for converting distance in m and azimuth angle (degrees, 
% 0 at north),
% 
% when given two locations in degrees latitude and longitude. 
% 
R = 6372795; % Earth radius in meters
lat1 = double(lat1) .* (pi/180);
lat2 = double(lat2) .* (pi/180);
lon1 = double(lon1) .* (pi/180);
lon2 = double(lon2) .* (pi/180);
d = acos( sin(lat1).*sin(lat2) + cos(lat1).*cos(lat2).*cos(lon2-lon1) ).*R;
mlat = (lat1+lat2)./2;
xd=sign(lon2-lon1).*R.*acos(sin(mlat).*sin(mlat)+cos(mlat).*cos(mlat).*cos(lon2-lon1));
yd=sign(lat2-lat1).*R.*acos(sin(lat2).*sin(lat1)+cos(lat2).*cos(lat1).*cos(lon1-lon1));



xd(~isreal(xd)) = NaN;
yd(~isreal(yd)) = NaN;

nnz(~isreal(xd(:)));
xd = real(xd);
yd = real(yd);

az = wrapTo360(atan2d(xd,yd));
