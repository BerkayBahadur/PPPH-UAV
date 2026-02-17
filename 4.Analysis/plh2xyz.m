function [x,y,z] = plh2xyz(lat,lon,h)

a  = 6378137.0;
f  = 1/298.257223563;
e2 = 2*f - f^2;

N = a./(sqrt(1 - e2*(sind(lat).^2)));

x = (N+h).*cosd(lat).*cosd(lon);
y = (N+h).*cosd(lat).*sind(lon);
z = ((1-e2)*N + h).*sind(lat);
end