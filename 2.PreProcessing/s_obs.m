function [data] = s_obs(data)

c = 299792458;

rn = size(data.rnx.obs.ep,1);
cn = size(data.rnx.obs.gps,2) + size(data.rnx.obs.glo,2) + size(data.rnx.obs.gal,2) + ...
    size(data.rnx.obs.bds,2);

pob = NaN(rn,cn);
sob = NaN(rn,cn,4);
sos.gps = cell(4,1); sos.glo = cell(4,1);
sos.gal = cell(4,1); sos.bds = cell(4,1);
fob = NaN(cn,1);

%GPS
prg1 = 'WCLX';
gp1 = NaN; gl1 = NaN;
for pi = 1:length(prg1)
    p1 = strcat('C1',prg1(pi));
    if any(strcmp(data.rnx.head.seq.gps,p1))
        gp1 = find(strcmp(data.rnx.head.seq.gps,p1));
        sos.gps{1} = p1;
        break
    end
end
for pi = 1:length(prg1)
    l1 = strcat('L1',prg1(pi));
    if any(strcmp(data.rnx.head.seq.gps,l1))
        gl1 = find(strcmp(data.rnx.head.seq.gps,l1));
        sos.gps{2} = l1;
        break
    end
end
prg2 = 'WSLX';
gp2 = NaN; gl2 = NaN;
for pi = 1:length(prg2)
    p2 = strcat('C2',prg2(pi));
    if any(strcmp(data.rnx.head.seq.gps,p2))
        gp2 = find(strcmp(data.rnx.head.seq.gps,p2));
        sos.gps{3} = p2;
        break
    end
end
for pi = 1:length(prg2)
    l2 = strcat('L2',prg2(pi));
    if any(strcmp(data.rnx.head.seq.gps,l2))
        gl2 = find(strcmp(data.rnx.head.seq.gps,l2));
        sos.gps{4} = l2;
        break
    end
end
if ~any(isnan([gp1,gl1,gp2,gl2]))
    [obb] = d_osbtype(data,'G',sos.gps);
    p1 = data.rnx.obs.gps(:,:,gp1) - (obb.ob1.*c)';
    l1 = data.rnx.obs.gps(:,:,gl1) - (obb.ob2.*c)';
    p2 = data.rnx.obs.gps(:,:,gp2) - (obb.ob3.*c)';
    l2 = data.rnx.obs.gps(:,:,gl2) - (obb.ob4.*c)';

    pob(:,1:32) = ~isnan(p1+l1+p2+l2);
    sob(:,1:32,1) = p1;
    sob(:,1:32,2) = l1;
    sob(:,1:32,3) = p2;
    sob(:,1:32,4) = l2;

    dod = obb.ob1 + obb.ob2 + obb.ob3 + obb.ob4;
    fob(1:32) = dod~=0;
end

% GLONASS
prr = 'PC';
rp1 = NaN; rl1 = NaN;
for pi = 1:length(prr)
    p1 = strcat('C1',prr(pi));
    if any(strcmp(data.rnx.head.seq.glo,p1))
        rp1 = find(strcmp(data.rnx.head.seq.glo,p1));
        sos.glo{1} = p1;
        break
    end
end
for pi = 1:length(prr)
    l1 = strcat('L1',prr(pi));
    if any(strcmp(data.rnx.head.seq.glo,l1))
        rl1 = find(strcmp(data.rnx.head.seq.glo,l1));
        sos.glo{2} = l1;
        break
    end
end
rp2 = NaN; rl2 = NaN;
for pi = 1:length(prr)
    p2 = strcat('C2',prr(pi));
    if any(strcmp(data.rnx.head.seq.glo,p2))
        rp2 = find(strcmp(data.rnx.head.seq.glo,p2));
        sos.glo{3} = p2;
        break
    end
end
for pi = 1:length(prr)
    l2 = strcat('L2',prr(pi));
    if any(strcmp(data.rnx.head.seq.glo,l2))
        rl2 = find(strcmp(data.rnx.head.seq.glo,l2));
        sos.glo{4} = l2;
        break
    end
end

if ~any(isnan([rp1,rl1,rp2,rl2]))
    [obb] = d_osbtype(data,'R',sos.glo);
    p1 = data.rnx.obs.glo(:,:,rp1) - (obb.ob1.*c)';
    l1 = data.rnx.obs.glo(:,:,rl1) - (obb.ob2.*c)';
    p2 = data.rnx.obs.glo(:,:,rp2) - (obb.ob3.*c)';
    l2 = data.rnx.obs.glo(:,:,rl2) - (obb.ob4.*c)';

    pob(:,33:58) = ~isnan(p1+l1+p2+l2);
    sob(:,33:58,1) = p1;
    sob(:,33:58,2) = l1;
    sob(:,33:58,3) = p2;
    sob(:,33:58,4) = l2;

    dod = obb.ob1 + obb.ob2 + obb.ob3 + obb.ob4;
    fob(33:58) = dod~=0;
end
% Galileo
pre1 = 'CX';
ep1 = NaN; el1 = NaN;
for pi = 1:length(pre1)
    p1 = strcat('C1',pre1(pi));
    if any(strcmp(data.rnx.head.seq.gal,p1))
        ep1 = find(strcmp(data.rnx.head.seq.gal,p1));
        sos.gal{1} = p1;
        break
    end
end
for pi = 1:length(pre1)
    l1 = strcat('L1',pre1(pi));
    if any(strcmp(data.rnx.head.seq.gal,l1))
        el1 = find(strcmp(data.rnx.head.seq.gal,l1));
        sos.gal{2} = l1;
        break
    end
end
ep2 = NaN; el2 = NaN;
pre5 = 'IQX';
for pi = 1:length(pre5)
    p2 = strcat('C5',pre5(pi));
    if any(strcmp(data.rnx.head.seq.gal,p2))
        ep2 = find(strcmp(data.rnx.head.seq.gal,p2));
        sos.gal{3} = p2;
        break
    end
end
for pi = 1:length(pre5)
    l2 = strcat('L5',pre5(pi));
    if any(strcmp(data.rnx.head.seq.gal,l2))
        el2 = find(strcmp(data.rnx.head.seq.gal,l2));
        sos.gal{4} = l2;
        break
    end
end

if ~any(isnan([ep1,el1,ep2,el2]))
    [obb] = d_osbtype(data,'E',sos.gal);
    p1 = data.rnx.obs.gal(:,:,ep1) - (obb.ob1.*c)';
    l1 = data.rnx.obs.gal(:,:,el1) - (obb.ob2.*c)';
    p2 = data.rnx.obs.gal(:,:,ep2) - (obb.ob3.*c)';
    l2 = data.rnx.obs.gal(:,:,el2) - (obb.ob4.*c)';

    pob(:,59:94) = ~isnan(p1+l1+p2+l2);
    sob(:,59:94,1) = p1;
    sob(:,59:94,2) = l1;
    sob(:,59:94,3) = p2;
    sob(:,59:94,4) = l2;

    dod = obb.ob1 + obb.ob2 + obb.ob3 + obb.ob4;
    fob(59:94) = dod~=0;
end

% BeiDou
prc = 'IQX';
cp1 = NaN; cl1 = NaN;
for pi = 1:length(prc)
    p1 = strcat('C2',prc(pi));
    if any(strcmp(data.rnx.head.seq.bds,p1))
        cp1 = find(strcmp(data.rnx.head.seq.bds,p1));
        sos.bds{1} = p1;
        break
    end
end
for pi = 1:length(prc)
    l1 = strcat('L2',prc(pi));
    if any(strcmp(data.rnx.head.seq.bds,l1))
        cl1 = find(strcmp(data.rnx.head.seq.bds,l1));
        sos.bds{2} = l1;
        break
    end
end
cp2 = NaN; cl2 = NaN;
for pi = 1:length(prc)
    p2 = strcat('C6',prc(pi));
    if any(strcmp(data.rnx.head.seq.bds,p2))
        cp2 = find(strcmp(data.rnx.head.seq.bds,p2));
        sos.bds{3} = p2;
        break
    end
end
for pi = 1:length(prc)
    l2 = strcat('L6',prc(pi));
    if any(strcmp(data.rnx.head.seq.bds,l2))
        cl2 = find(strcmp(data.rnx.head.seq.bds,l2));
        sos.bds{4} = l2;
        break
    end
end
if ~any(isnan([cp1,cl1,cp2,cl2]))
    [obb] = d_osbtype(data,'C',sos.bds);
    p1 = data.rnx.obs.bds(:,:,cp1) - (obb.ob1.*c)';
    l1 = data.rnx.obs.bds(:,:,cl1) - (obb.ob2.*c)';
    p2 = data.rnx.obs.bds(:,:,cp2) - (obb.ob3.*c)';
    l2 = data.rnx.obs.bds(:,:,cl2) - (obb.ob4.*c)';

    pob(:,95:156) = ~isnan(p1+l1+p2+l2);
    sob(:,95:156,1) = p1;
    sob(:,95:156,2) = l1;
    sob(:,95:156,3) = p2;
    sob(:,95:156,4) = l2;

    dod = obb.ob1 + obb.ob2 + obb.ob3 + obb.ob4;
    fob(95:156) = dod~=0;
end
%
data.rnx.pob = pob;
data.rnx.sob = sob;
data.rnx.sos = sos;
data.rnx.fob = fob;

end