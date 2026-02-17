function [sapc] = c_satapc(s_xyz,r_xyz,sun_xyz,sapc,sno,opt)

l   = s_xyz - r_xyz;
los = l./norm(l);

k = (-1).*(s_xyz./(norm(s_xyz)));
rs  = sun_xyz - s_xyz;
e   = rs./(norm(rs));
j   = cross(k,e)./norm(cross(k,e));
i   = cross(j,k);
sf  = [i; j; k];

neu = sapc.neu(sno,:,opt)';
rk = sf\neu;

l2 = r_xyz - s_xyz; 
ls2 = l2./norm(l2);
zen = 90 - asind(dot(ls2,k));
az = mod(atan2d(dot(ls2,i),dot(ls2,j)),360);

if ~isempty(sapc.adp{sno,opt})
    arr = sapc.adp{sno,opt};
    rs = find(az-arr(:,1)>0,1,'last'); 
    rf = rs + 1;
    xs = sapc.daz(sno,2):sapc.daz(sno,4):sapc.daz(sno,3);
    ys = mean(arr(rs:rf,2:end));
    pcv = interp1(xs,ys,zen);
else
    arr = sapc.nad{sno,opt};
    xs = sapc.daz(sno,2):sapc.daz(sno,4):sapc.daz(sno,3);
    pcv = interp1(xs,arr,zen);
end

sapc = (dot(rk,los) + pcv)./1000;

end