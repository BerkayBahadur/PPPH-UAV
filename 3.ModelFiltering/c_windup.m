function [wup] = c_windup(rec,sat,sun,prev)

% unit vector from satellite to sun
sat2sun = (sun - sat)/norm(sun - sat,'fro');
% axes in satellite body-fixed system
ez = -sat/norm(sat,'fro');
ey0 = cross(ez,sat2sun);
ey = ey0/norm(ey0,'fro');
ex = cross(ey,ez);

% from XYZ to phi, lambda, h
[elip] = xyz2plh(rec,0);
phi = elip(1);
lam = elip(2);

% north in local
xr = [-sin(phi)*cos(lam) -sin(phi)*sin(lam) cos(phi)];
% east in local
yr = [-sin(lam) cos(lam) 0];

% line of sight (unit vector) from satellite to receiver
p0 = rec - sat;
p = p0./norm(p0);

% receiver east and north directions
ar = yr; br = xr;
% satellite i and j directions
as = ex; bs = ey;

% dipoles for receiver and satellite
Dr = ar - p*dot(p,ar) + cross(p,br);
Ds = as - p*dot(p,as) - cross(p,bs);

cosphi = dot(Ds,Dr)/(norm(Ds)*norm(Dr));
if cosphi>1, cosphi=1; end

y = dot(p,cross(Ds,Dr));
dphi = sign(y)*acos(cosphi);

if prev == 0
    N = 0;
else
    dphi_prev = prev*2*pi; % radian
    N = round((dphi_prev - dphi)/(2*pi));
end

Dphi = dphi + N*2*pi;

wup = Dphi/(2*pi); %in cycle

end