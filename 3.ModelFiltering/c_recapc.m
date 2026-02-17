function [rapc] = c_recapc(s_xyz,r_xyz,rapc,kn,az,elev)

l   = r_xyz - s_xyz;
los = l./norm(l);

[elip] = xyz2plh(r_xyz,0);
lat = elip(1);
lon = elip(2);

ori = [(-sin(lon)) (-cos(lon)*sin(lat)) (cos(lon)*cos(lat));...
        (cos(lon)) (-sin(lon)*sin(lat)) (sin(lon)*cos(lat));...
               (0) (cos(lat))           (sin(lat))];
           
f = rapc.neu(:,:,kn);
if (sum(f==0)==3)
    if rem(kn,2)==1
        f = rapc.neu(:,:,1);
    else
        f = rapc.neu(:,:,2);
    end
end

c = [f(2);
     f(1);
     f(3)];

p = ori*c;
if size(los,1)~=size(p,1)
    los = los';
end

zen = 90 - elev;
if ~isempty(rapc.adp{kn})
    arr = rapc.adp{kn};
else
    if rem(kn,2)==1
        arr = rapc.adp{1};
    else
        arr = rapc.adp{2};
    end
end

if ~isempty(arr)
    rs = find(az-arr(:,1)>0,1,'last'); rf = rs + 1;
    xs = rapc.daz(2):rapc.daz(4):rapc.daz(3);
    ys = mean(arr(rs:rf,2:end));
    pcv = interp1(xs,ys,zen);
else
    pcv = 0;
end


rapc = (dot(p,los) + pcv)./1000;

end