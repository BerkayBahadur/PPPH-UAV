function [xyz] = NEU2xyz(DnDeDu,lat,lon)


R = [-cos(lon)*sin(lat) -sin(lon)*sin(lat) cos(lat);
     -sin(lon)          cos(lon)           0       ;
     cos(lon)*cos(lat)  sin(lon)*cos(lat)  sin(lat)];
xyz = ((R^-1)*DnDeDu')';