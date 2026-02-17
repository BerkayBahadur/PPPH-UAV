function [data] = c_satpos(data)

c = 299792458;
we=7.2921151467e-5;

sp3int = data.sp3.sp3int;
clkint = data.clk.clkint;
rec    = data.rnx.head.rec.pos;
if size(rec,1) ~= 1
    rec = rec';
end

en = size(data.rnx.pob,1);
sn = size(data.rnx.pob,2);
crt = NaN(en,7,sn);
tfs = NaN(en,sn);

for i=1:en
    ls = find(data.rnx.pob(i,:) == 1);
    for k=ls
        if k<33
            kn = k;
            cs = data.clk.gps(:,kn);
            xs = data.sp3.gps(:,1,kn);
            ys = data.sp3.gps(:,2,kn);
            zs = data.sp3.gps(:,3,kn);

        elseif k<59
            kn = k - 32;
            cs = data.clk.glo(:,kn);
            xs = data.sp3.glo(:,1,kn);
            ys = data.sp3.glo(:,2,kn);
            zs = data.sp3.glo(:,3,kn);

        elseif k<95
            kn = k - 58;
            cs = data.clk.gal(:,kn);
            xs = data.sp3.gal(:,1,kn);
            ys = data.sp3.gal(:,2,kn);
            zs = data.sp3.gal(:,3,kn);

        else
            kn = k - 94;
            cs = data.clk.bds(:,kn);
            xs = data.sp3.bds(:,1,kn);
            ys = data.sp3.bds(:,2,kn);
            zs = data.sp3.bds(:,3,kn);
            
        end

        tof = data.rnx.sob(i,k,1)/c;
        nep = data.rnx.obs.ep(i,1) - tof;
        dt  = entrp(nep,clkint,cs);
        if ~isnan(dt)
            dt  = entrp(nep - dt,clkint,cs);
        else
            data.rnx.pob(i,k) = 0;
            continue
        end
        nep = nep - dt;
        tfs(i,k) = tof + dt;

        [x,vx] = entrp(nep,sp3int,xs);
        [y,vy] = entrp(nep,sp3int,ys);
        [z,vz] = entrp(nep,sp3int,zs);
        [dt,~] = entrp(nep,clkint,cs);

        R  = [x y z];
        tf     = norm(R - rec)./c;
        er_ang = rad2deg(tf*we);
        pos    = rotation(R,er_ang,3);

        if any(isnan(pos))
            data.rnx.pob(i,k) = 0;
        else
            crt(i,1,k) = pos(1);
            crt(i,2,k) = pos(2);
            crt(i,3,k) = pos(3);
            crt(i,4,k) = vx;
            crt(i,5,k) = vy;
            crt(i,6,k) = vz;
            crt(i,7,k) = dt;
        end
    end
end

data.sat.crt = crt;
data.sat.tfs = tfs;
end