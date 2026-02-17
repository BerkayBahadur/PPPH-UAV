function [data] = c_osb(data)

sn = size(data.rnx.pob,2);
c   = 299792458;% m/s
[frqs,~] = frequencies;
tosb = NaN(sn,3);

for k=1:sn
    f1 = frqs(k,1); f2 = frqs(k,2);
    wl1 =  f1/(f1-f2); wl2 = -f2/(f1-f2);
    nl1 =  f1/(f1+f2); nl2 =  f2/(f1+f2);
    if1 = (f1^2)/(f1^2-f2^2); if2 = -(f2^2)/(f1^2-f2^2);
    if k<33
        prn = k;
        tosb(k,1) = (wl1*data.osb.gps.l1w(prn) + wl2*data.osb.gps.l2w(prn))*c; % WL phase (m)
        tosb(k,2) = (nl1*data.osb.gps.c1w(prn) + nl2*data.osb.gps.c2w(prn))*c; % NL code (m)
        tosb(k,3) = (if1*data.osb.gps.l1w(prn) + if2*data.osb.gps.l2w(prn))*c; % NL phase (m)

    elseif k<59
        prn = k - 32;
        tosb(k,1) = (wl1*data.osb.glo.l1p(prn) + wl2*data.osb.glo.l2p(prn))*c;
        tosb(k,2) = (nl1*data.osb.glo.c1p(prn) + nl2*data.osb.glo.c2p(prn))*c;
        tosb(k,3) = (if1*data.osb.glo.l1p(prn) + if2*data.osb.glo.l2p(prn))*c;

    elseif k<95
        prn = k - 58;
        tosb(k,1) = (wl1*data.osb.gal.l1c(prn) + wl2*data.osb.gal.l5q(prn))*c;
        tosb(k,2) = (nl1*data.osb.gal.c1c(prn) + nl2*data.osb.gal.c5q(prn))*c;
        tosb(k,3) = (if1*data.osb.gal.l1c(prn) + if2*data.osb.gal.l5q(prn))*c;

    else
        prn = k - 94;
        tosb(k,1) = (wl1*data.osb.bds.l2i(prn) + wl2*data.osb.bds.l6i(prn))*c;
        tosb(k,2) = (nl1*data.osb.bds.c2i(prn) + nl2*data.osb.bds.c6i(prn))*c;
        tosb(k,3) = (if1*data.osb.bds.l2i(prn) + if2*data.osb.bds.l6i(prn))*c;

    end
end

data.osb.tosb = tosb;

end