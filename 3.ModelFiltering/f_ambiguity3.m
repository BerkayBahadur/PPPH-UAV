function [refs,stes,sdwl,fxwl,fxll,fpos] = f_ambiguity3(refs,meas,data,stes,epn,sdwl,bp,sls,amb,pfxll)

c   = 299792458;% m/s
[freq,~] = frequencies;

fxwl = NaN(1,156);
fxll = NaN(1,156);
%%
als = unique(meas(:,4));
gls = als(als<33);
rls = als(als>32 & als<59);
els = als(als>58 & als<95);
c2ls = als(als>94 & als<111);
c3ls = als(als>110);

ca = 30;
if isnan(refs(1)) || sum(meas(:,4)==refs(1))==0 ||...
        (sum(meas(:,4)==refs(1))>0 && unique(meas(meas(:,4)==refs(1),26))<ca)
    elv = NaN(size(gls,1),1);
    for k=1:size(gls,1)
        elv(k) = unique(meas(meas(:,4)==gls(k),26));
    end
    if gls(elv==max(elv))~=refs(1)
        refs(1) = gls(elv==max(elv));
        stes(1) = epn;
    end
end

if isnan(refs(2)) || sum(meas(:,4)==refs(2))==0 ||...
        (sum(meas(:,4)==refs(2))>0 && unique(meas(meas(:,4)==refs(2),26))<ca)
    elv = NaN(size(rls,1),1);
    for k=1:size(rls,1)
        elv(k) = unique(meas(meas(:,4)==rls(k),26));
    end
    if rls(elv==max(elv))~=refs(2)
        refs(2) = rls(elv==max(elv));
        stes(2) = epn;
    end
end

if isnan(refs(3)) || sum(meas(:,4)==refs(3))==0 ||...
        (sum(meas(:,4)==refs(3))>0 && unique(meas(meas(:,4)==refs(3),26))<ca)
    elv = NaN(size(els,1),1);
    for k=1:size(els,1)
        elv(k) = unique(meas(meas(:,4)==els(k),26));
    end
    if els(elv==max(elv))~=refs(3)
        refs(3) = els(elv==max(elv));
        stes(3) = epn;
    end
end

if isnan(refs(4)) || sum(meas(:,4)==refs(4))==0 ||...
        (sum(meas(:,4)==refs(4))>0 && unique(meas(meas(:,4)==refs(4),26))<ca)
    elv = NaN(size(c2ls,1),1);
    for k=1:size(c2ls,1)
        elv(k) = unique(meas(meas(:,4)==c2ls(k),26));
    end
    if c2ls(elv==max(elv))~=refs(4)
        refs(4) = c2ls(elv==max(elv));
        stes(4) = epn;
    end
end

if isnan(refs(5)) || sum(meas(:,4)==refs(5))==0 ||...
        (sum(meas(:,4)==refs(5))>0 && unique(meas(meas(:,4)==refs(5),26))<ca)
    elv = NaN(size(c3ls,1),1);
    for k=1:size(c3ls,1)
        elv(k) = unique(meas(meas(:,4)==c3ls(k),26));
    end
    if c3ls(elv==max(elv))~=refs(5)
        refs(5) = c3ls(elv==max(elv));
        stes(5) = epn;
    end
end
%%
hmw = NaN(1,156);
for i=1:length(als)
    whl = meas(meas(:,4)==als(i),:);
    f1 = freq(als(i),1); f2 = freq(als(i),2);
    p1 = whl(1,6) - whl(1,16); p2 = whl(2,6) - whl(2,16);
    l1 = whl(3,6) - whl(3,16); l2 = whl(4,6) - whl(4,16);

    lwl = (l1*f1 - l2*f2)/(f1-f2);
    pnl = (p1*f1 + p2*f2)/(f1+f2);
    lamwl = c/(f1-f2);
    hmw(als(i)) = (lwl - pnl)/(lamwl);
end

sd_hmw = NaN(1,156);
if ~isempty(gls)
    sd_hmw(gls)  = hmw(refs(1)) - hmw(gls);
end
if ~isempty(rls)
    sd_hmw(rls)  = hmw(refs(2)) - hmw(rls);
end
if ~isempty(els)
    sd_hmw(els)  = hmw(refs(3)) - hmw(els);
end
if ~isempty(c2ls)
    sd_hmw(c2ls) = hmw(refs(4)) - hmw(c2ls);
end
if ~isempty(c3ls)
    sd_hmw(c3ls) = hmw(refs(5)) - hmw(c3ls);
end
sdwl(epn,:) = sd_hmw;

%
cen0 = (120/data.rnx.head.time.int) - 1;
cen1 = (240/data.rnx.head.time.int) - 1;
% 
cfc = 0.25;
for i=1:length(als)
    snn = als(i);
    if snn<33
        rs = stes(1);
    elseif snn<59
        rs = stes(2);
    elseif snn<95
        rs = stes(3);
    elseif snn<111
        rs = stes(4);
    else
        rs = stes(5);
    end
    % 
    if any(isnan(sdwl(1:epn,snn)))
        as = find(isnan(sdwl(1:epn,snn)),1,"last") + 1;
    else
        as = 1;
    end
    % 
    aln = length(as:epn);
    if aln>cen0
        % 
        if as<rs
            ds = sdwl(rs:epn,snn);
        else
            ds = sdwl(as:epn,snn);
        end
    
        if length(ds)>cen1
            fwl = mean(ds(end-cen1:end),'omitmissing');
        else
            fwl = mean(ds,'omitmissing');
        end

        % 
        if abs(fwl - round(fwl))<cfc
            fxwl(snn) = round(fwl);
        end
    end
end

% 
lls = find(~isnan(fxwl) & fxwl~=0);
lls(data.rnx.fob(lls)==0) = [];
nlsd = NaN(length(lls),1); nsds = NaN(length(lls),1);
xxk = amb.xk(bp+1:end);
for i=1:length(lls)
    sno = lls(i);
    f1 = freq(sno,1); f2 = freq(sno,2);
    w1 = c/f1;
    ins = find(sls==sno);

    if sno<33
        inr = sls==refs(1);
        nsd = xxk(inr) - xxk(ins); 
        nsds(i) = nsd;
        nlsd(i) = nsd*((f1+f2)/(w1*f1)) - fxwl(sno)*(f2/(f1-f2));

    elseif sno<59
        inr = sls==refs(2);
        nsd = xxk(inr) - xxk(ins); 
        nsds(i) = nsd;
        nlsd(i) = nsd*((f1+f2)/(w1*f1)) - fxwl(sno)*(f2/(f1-f2));

    elseif sno<95
        inr = sls==refs(3);
        nsd = xxk(inr) - xxk(ins); 
        nsds(i) = nsd;
        nlsd(i) = nsd*((f1+f2)/(w1*f1)) - fxwl(sno)*(f2/(f1-f2));

    elseif sno<111
        inr = sls==refs(4);
        nsd = xxk(inr) - xxk(ins); 
        nsds(i) = nsd;
        nlsd(i) = nsd*((f1+f2)/(w1*f1)) - fxwl(sno)*(f2/(f1-f2));

    else
        inr = sls==refs(5);
        nsd = xxk(inr) - xxk(ins); 
        nsds(i) = nsd;
        nlsd(i) = nsd*((f1+f2)/(w1*f1)) - fxwl(sno)*(f2/(f1-f2));

    end
end

Qnn = amb.pk(bp+1:end,bp+1:end);
nsls = sls;
if ~any(lls<33)
    del = nsls<33;
    nsls(del) = [];
    Qnn(del,:) = [];
    Qnn(:,del) = [];
end
if ~any(lls>32 & lls<59)
    del = nsls>32 & nsls<59;
    nsls(del) = [];
    Qnn(del,:) = [];
    Qnn(:,del) = [];
end
if ~any(lls>58 & lls<95)
    del = nsls>58 & nsls<95;
    nsls(del) = [];
    Qnn(del,:) = [];
    Qnn(:,del) = [];
end
if ~any(lls>94)
    del = nsls>94;
    nsls(del) = [];
    Qnn(del,:) = [];
    Qnn(:,del) = [];
end

ck = -eye(length(nsls));
for k=1:length(refs)
    ins = find(nsls==refs(k));
    switch k
        case 1
            ll = nsls<33;
            ck(ll,ins) = 1;
            ck(ins,ins) = 0;
        case 2
            ll = nsls>32 & nsls<59;
            ck(ll,ins) = 1;
            ck(ins,ins) = 0;
        case 3
            ll = nsls>58 & nsls<95;
            ck(ll,ins) = 1;
            ck(ins,ins) = 0;
        case 4
            ll = nsls>94 & nsls<111;
            ck(ll,ins) = 1;
            ck(ins,ins) = 0;
        case 5
            ll = nsls>110;
            ck(ll,ins) = 1;
            ck(ins,ins) = 0;
    end

end

Qn = ck*Qnn*ck';
rem = ~ismember(nsls,lls);
Qn(rem,:) = [];
Qn(:,rem) = [];

if length(nlsd)>6
    try
        [fxnl,~,~,~,~,~] = LAMBDA(nlsd,Qn,5,'P0',0.995);
    catch
        fxnl = LAMBDA(nlsd,Qn, 1);
    end
    
    bs = fxnl(:,1) - floor(fxnl(:,1)) == 0;
    bfxnl = fxnl(bs,1);
    bls = lls(bs);
    %
    ddel = zeros(length(bfxnl),1);
    if any(~isnan(pfxll))
        for kk=1:length(bfxnl)
            snk = bls(kk);
            if snk<33
                ats = stes(1);
            elseif snk<59
                ats = stes(2);
            elseif snk<95
                ats = stes(3);
            elseif snk<111
                ats = stes(4);
            else
                ats = stes(5);
            end
            if ~isnan(pfxll(snk)) && ~isequal(bfxnl(kk),pfxll(snk)) && ~isequal(epn,ats)
                ddel(kk) = 1;
            end
        end
        %
        if any(ddel==1)
            bfxnl(ddel==1) = [];
            bls(ddel==1) = [];
        end
    end
    %
    nif = NaN(1,length(bfxnl));
    for i=1:length(nif)
        prn = bls(i);
        f1 = freq(prn,1); f2 = freq(prn,2);
        w1 = c/f1;
    
        nif(i) = ((w1*f1*f2)/(f1^2-f2^2))*fxwl(prn) + ((w1*f1)/(f1+f2))*(bfxnl(i,1));
        fxll(prn) = bfxnl(i,1);
    end
    
    Lk = [amb.vk;nif'];
    
    A = amb.hk;
    A(:,4:bp) = [];
    
    Ad = zeros(length(bls),size(A,2));
    for i=1:length(bls)
        prn = bls(i);
        if prn<33
            inp = find(sls==prn);
            inr = find(sls==refs(1));
            Ad(i,3+inp) = -1; Ad(i,3+inr) = 1;
        elseif prn<59
            inp = find(sls==prn);
            inr = find(sls==refs(2));
            Ad(i,3+inp) = -1; Ad(i,3+inr) = 1;
        elseif prn<95
            inp = find(sls==prn);
            inr = find(sls==refs(3));
            Ad(i,3+inp) = -1; Ad(i,3+inr) = 1;
        elseif prn<111
            inp = find(sls==prn);
            inr = find(sls==refs(4));
            Ad(i,3+inp) = -1; Ad(i,3+inr) = 1;
        else
            inp = find(sls==prn);
            inr = find(sls==refs(5));
            Ad(i,3+inp) = -1; Ad(i,3+inr) = 1;
        end
    
    end
    
    Ak = [A;Ad];
    
    Pfl = inv(amb.rk);
    
    
    glls = bls(bls<33);
    if ~isempty(glls)
        [Pg] = c_weight(amb.rk,sls,glls,refs(1));
    else
        Pg = [];
    end
    
    rlls = bls(bls>32 & bls<59);
    if ~isempty(rlls)
        [Pr] = c_weight(amb.rk,sls,rlls,refs(2));
    else
        Pr = [];
    end
    
    ells = bls(bls>58 & bls<95);
    if ~isempty(ells)
        [Pe] = c_weight(amb.rk,sls,ells,refs(3));
    else
        Pe = [];
    end
    
    c2lls = bls(bls>94 & bls<111);
    if ~isempty(c2lls)
        [Pc2] = c_weight(amb.rk,sls,c2lls,refs(4));
    else
        Pc2 = [];
    end
    
    c3lls = bls(bls>110);
    if ~isempty(c3lls)
        [Pc3] = c_weight(amb.rk,sls,c3lls,refs(5));
    else
        Pc3 = [];
    end
    
    Pk = blkdiag(Pfl,Pg,Pr,Pe,Pc2,Pc3);
    
    N = Ak'*Pk*Ak;              % Normal Equation Matrix
    Qxx = pinv(N);              % Cofactor Matrix of Parameters
    dx = Qxx*Ak'*Pk*Lk;   	    % Adjusted Parameters
    
    fpos = amb.x1(1:3) + dx(1:3);
else
    fpos = amb.xk(1:3);
end

end

function [Pg] = c_weight(Rk,sls,slls,refn)

nn = length(slls);
Qt = eye(nn);
for i=1:length(slls)
    il = find(sls==slls(i));
    Qt(i,i) = Rk(2*il,2*il);
end
ir = find(sls==refn);
Qt(end+1,end+1) = Rk(2*ir,2*ir);
C = zeros(nn,nn+1);
C(:,end) = 1;       % Position of reference-satellite
C(:,1:end-1) = -eye(nn);
Pg = (C*Qt*C')^(-1)*1000^2;      % 1000 seems to be a good choice

end
