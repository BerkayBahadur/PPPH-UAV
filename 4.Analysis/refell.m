function [a,b,e2,finv]=refell(type)


type=upper(type);
if (type=='CLK66' | type=='NAD27')
  a=6378206.4;
  finv=294.9786982;
elseif type=='GRS67'
  a=6378160.0;
  finv=298.247167427;
elseif (type=='GRS80' | type=='NAD83')
  a=6378137.0;
  finv=298.257222101;
elseif (type=='WGS72')
  a=6378135.0;
  finv=298.26;
elseif (type=='WGS84')
  a=6378137.0;
  finv=298.257223563;
elseif type=='ATS77'
  a=6378135.0;
  finv=298.257;
elseif type=='KRASS'
  a=6378245.0;
  finv=298.3;
elseif type=='INTER'
  a=6378388.0;
  finv=297.0;
elseif type=='MAIRY'
  a=6377340.189;
  finv=299.3249646;
elseif type=='TOPEX'
  a=6378136.3;
  finv=298.257;
end
f=1/finv;
b=a*(1-f);
e2=1-(1-f)^2;
