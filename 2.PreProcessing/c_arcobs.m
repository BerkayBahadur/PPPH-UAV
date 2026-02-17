function [data] = c_arcobs(data)

[arc] = d_arcs(data);
en = size(data.rnx.pob,1);
sn = size(data.rnx.pob,2);
pob = zeros(en,sn);

for k=1:sn
    karc = arc{k};
    for s=1:size(karc,1)
        st = karc(s,1); fn = karc(s,2);
        pob(st:fn,k) = data.rnx.pob(st:fn,k);
    end
end

data.rnx.pob = pob;

end