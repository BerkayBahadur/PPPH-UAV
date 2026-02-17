function [Trop,Mwet,Mn,Me,ZHD,ZWD_ap] = c_trop(rec,sat,dmjd,p,ah,aw,Tm,e,la)

% Ellipsoidal Geodetic Coordinates
[ellp] = xyz2plh(rec,0);
dlat = ellp(1); % radian
dlon = ellp(2); % radian
hell = ellp(3); % meter

% Local azimuth and elevation angles
[Az,Elv] = local(rec,sat,0);

% Saastamoinen Model for dry tropospheric part
f = 0.0022768*p;
k = 1 - (0.00266*cos(2*dlat)) - ((0.28*10^-6)*hell);
ZHD = f/k;

% Saastamoinen Model for wet tropospheric part
Rd = 287.0464; % J/(K*kg) 1J = kg*m2/s2
gm = 9.80665;  % m/s2
k1 = 77.6890;  % K/hPa
k2 = 71.2952;  % K/hPa
k3 = 375463;   % K2/hPa
% These coefficents are taken from Path Delays in the Neutral Atmosphere,
% T. Nilsson, J. Böhm, et. al, 2013.
mw = 18.0152; % g/mol molar massof water
md = 28.9644; % g/mol molar mass of dry air
k2d = k2 - (k1*(mw/md)); % K/hPa

ZWD_ap = (10^-6)*(k2d + (k3/Tm))*((Rd*e)/(gm*(la+1)));% meter

[mfh,mfw] = c_vmf3ht(ah,aw,dmjd,dlat,dlon,hell,(pi/2 - Elv));

Trop = mfh*ZHD + mfw*ZWD_ap;
Mwet = mfw;

Mg = 1/((tan(Elv)*sin(Elv))+0.0032);
Mn = Mg*cos(Az);
Me = Mg*sin(Az);

end