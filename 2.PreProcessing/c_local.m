function [data] = c_local(data,elvang)

en = size(data.rnx.pob,1);
sn = size(data.rnx.pob,2);

rec = data.rnx.head.rec.pos;

lcl = NaN(en,sn,2);

for i=1:en
    ls = find(data.rnx.pob(i,:) == 1);
    for k=ls
        if data.rnx.pob(i,k) == 1
            sat = data.sat.crt(i,1:3,k);
            [az,el] = local(rec,sat,1);
            lcl(i,k,1) = az;
            lcl(i,k,2) = el;
            if el<elvang
                data.rnx.pob(i,k) = 0;
            end
        end
    end
end

data.sat.lcl = lcl;

end