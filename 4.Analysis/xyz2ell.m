function [lat,lon,h]=xyz2ell(X,Y,Z,a,e2)


if nargin ~= 3 & nargin ~= 5
  warning('Incorrect number of input arguments');
  return
end
if nargin == 3
  [a,b,e2]=refell('grs80');
end

elat=1.e-12;
eht=1.e-5;

p=sqrt(X.*X+Y.*Y);
lat=atan2(Z,p.*(1-e2));
h=0;
dh=1;
dlat=1;

while sum(dlat>elat) | sum(dh>eht)
  lat0=lat;
  h0=h;
  v=a./sqrt(1-e2.*sin(lat).*sin(lat));
  h=p./cos(lat)-v;
  lat=atan2(Z, p.*(1-e2.*v./(v+h)));
  dlat=abs(lat-lat0);
  dh=abs(h-h0);
end
lon=atan2(Y,X);
