function [ xk,pk,kof,res,amb ] = f_rkalman(sno,sls,pno,meas,x1k,p1k,sit,bp,freq,options)

mno = 2*sno;
Hk = zeros(mno,pno);
Zk = zeros(mno,  1);
Ck = zeros(mno,  1);
Rk = eye(mno);
Rw = zeros(mno,  1);

for k=1:sno
    sat = meas((4*k)-3,8:10);
    rho = norm(sat - x1k(1:3,1)');
    
    Jx = (x1k(1,1) - sat(1,1))/rho;
    Jy = (x1k(2,1) - sat(1,2))/rho;
    Jz = (x1k(3,1) - sat(1,3))/rho;
    Jw = meas((4*k - 3),28);
    Jt = 1;
    Jn = 1;

    s = (2*k - 1);
    f = (2*k);
    Hk(s:f,1     ) = Jx;
    Hk(s:f,2     ) = Jy;
    Hk(s:f,3     ) = Jz;
    Hk(s:f,4     ) = Jw;
    Hk(f  ,(bp+k)) = Jn;

    if sit==1 || sit==2 || sit==3 % G R E
        Hk(s:f,5) = Jt;
        
    elseif sit==4 % C
        if sls(k)<111
            Hk(s:f,5) = Jt;
        else
            Hk(s:f,6) = Jt;
        end
        
    elseif sit==5 || sit==6 % GR GE
        if sls(k)<33
            Hk(s:f,5) = Jt;
        else
            Hk(s:f,6) = Jt;
        end

    elseif sit==7 % GC
        if sls(k)<33
            Hk(s:f,5) = Jt;
        elseif sls(k)<111
            Hk(s:f,6) = Jt;
        else
            Hk(s:f,7) = Jt;
        end

    elseif sit==8 % RE
        if sls(k)<59
            Hk(s:f,5) = Jt;
        else
            Hk(s:f,6) = Jt;
        end

    elseif sit==9 || sit==10 %RC EC
        if sls(k)<95
            Hk(s:f,5) = Jt;
        elseif sls(k)<111
            Hk(s:f,6) = Jt;
        else
            Hk(s:f,7) = Jt;
        end

    elseif sit==11 %GRE
        if sls(k)<33
            Hk(s:f,5) = Jt;
        elseif sls(k)<59
            Hk(s:f,6) = Jt;
        else
            Hk(s:f,7) = Jt;
        end

    elseif sit==12 %GRC
        if sls(k)<33
            Hk(s:f,5) = Jt;
        elseif sls(k)<59
            Hk(s:f,6) = Jt;
        elseif sls(k)<111
            Hk(s:f,7) = Jt;
        else
            Hk(s:f,8) = Jt;
        end
     elseif sit==13 %GEC
        if sls(k)<33
            Hk(s:f,5) = Jt;
        elseif sls(k)<95
            Hk(s:f,6) = Jt;
        elseif sls(k)<111
            Hk(s:f,7) = Jt;
        else
            Hk(s:f,8) = Jt;
        end

    elseif sit==14 %REC
        if sls(k)<59
            Hk(s:f,5) = Jt;
        elseif sls(k)<95
            Hk(s:f,6) = Jt;
        elseif sls(k)<111
            Hk(s:f,7) = Jt;
        else
            Hk(s:f,8) = Jt;
        end
    elseif sit==15 %GREC
        if sls(k)<33
            Hk(s:f,5) = Jt;
        elseif sls(k)<59
            Hk(s:f,6) = Jt;
        elseif sls(k)<95
            Hk(s:f,7) = Jt;
        elseif sls(k)<111
            Hk(s:f,8) = Jt;
        else
            Hk(s:f,9) = Jt;
        end
    end

    % iono-free measurement
    f1 = freq(sls(k),1); f2 = freq(sls(k),2);
    c1 = (f1^2)/(f1^2-f2^2); c2 = (f2^2)/(f1^2-f2^2);
    Zk(s,1) = c1*meas((4*k - 3),6) - c2*meas((4*k - 2),6);
    Zk(f,1) = c1*meas((4*k - 1),6) - c2*meas((4*k    ),6);
    % iono-free corrections for code and phase observations
    p1c = (meas((4*k - 3),7));
    p2c = (meas((4*k - 2),7));
    l1c = (meas((4*k - 1),7));
    l2c = (meas((4*k    ),7));
    pc = c1*p1c - c2*p2c;
    lc = c1*l1c - c2*l2c;

    if sit==1 || sit==2 || sit==3 % G R E
        Ck(s,1) = rho + pc + (Jw*x1k(4)) + (Jt*x1k(5));
        Ck(f,1) = rho + lc + (Jw*x1k(4)) + (Jt*x1k(5)) + (Jn*x1k(bp+k));

    elseif sit==4 % C
        if sls(k)<111
            Ck(s,1) = rho + pc + (Jw*x1k(4)) + (Jt*x1k(5));
            Ck(f,1) = rho + lc + (Jw*x1k(4)) + (Jt*x1k(5)) + (Jn*x1k(bp+k));
        else
            Ck(s,1) = rho + pc + (Jw*x1k(4)) + (Jt*x1k(6));
            Ck(f,1) = rho + lc + (Jw*x1k(4)) + (Jt*x1k(6)) + (Jn*x1k(bp+k));
        end
        
    elseif sit==5 || sit==6 % GR GE
        if sls(k)<33
            Ck(s,1) = rho + pc + (Jw*x1k(4)) + (Jt*x1k(5));
            Ck(f,1) = rho + lc + (Jw*x1k(4)) + (Jt*x1k(5)) + (Jn*x1k(bp+k));
        else
            Ck(s,1) = rho + pc + (Jw*x1k(4)) + (Jt*x1k(6));
            Ck(f,1) = rho + lc + (Jw*x1k(4)) + (Jt*x1k(6)) + (Jn*x1k(bp+k));
        end

    elseif sit==7 % GC
        if sls(k)<33
            Ck(s,1) = rho + pc + (Jw*x1k(4)) + (Jt*x1k(5));
            Ck(f,1) = rho + lc + (Jw*x1k(4)) + (Jt*x1k(5)) + (Jn*x1k(bp+k));
        elseif sls(k)<111
            Ck(s,1) = rho + pc + (Jw*x1k(4)) + (Jt*x1k(6));
            Ck(f,1) = rho + lc + (Jw*x1k(4)) + (Jt*x1k(6)) + (Jn*x1k(bp+k));
        else
            Ck(s,1) = rho + pc + (Jw*x1k(4)) + (Jt*x1k(7));
            Ck(f,1) = rho + lc + (Jw*x1k(4)) + (Jt*x1k(7)) + (Jn*x1k(bp+k));
        end

    elseif sit==8 % RE
        if sls(k)<59
            Ck(s,1) = rho + pc + (Jw*x1k(4)) + (Jt*x1k(5));
            Ck(f,1) = rho + lc + (Jw*x1k(4)) + (Jt*x1k(5)) + (Jn*x1k(bp+k));
        else
            Ck(s,1) = rho + pc + (Jw*x1k(4)) + (Jt*x1k(6));
            Ck(f,1) = rho + lc + (Jw*x1k(4)) + (Jt*x1k(6)) + (Jn*x1k(bp+k));
        end

    elseif sit==9 || sit==10 %RC EC
        if sls(k)<95
            Ck(s,1) = rho + pc + (Jw*x1k(4)) + (Jt*x1k(5));
            Ck(f,1) = rho + lc + (Jw*x1k(4)) + (Jt*x1k(5)) + (Jn*x1k(bp+k));
        elseif sls(k)<111
            Ck(s,1) = rho + pc + (Jw*x1k(4)) + (Jt*x1k(6));
            Ck(f,1) = rho + lc + (Jw*x1k(4)) + (Jt*x1k(6)) + (Jn*x1k(bp+k));
        else
            Ck(s,1) = rho + pc + (Jw*x1k(4)) + (Jt*x1k(7));
            Ck(f,1) = rho + lc + (Jw*x1k(4)) + (Jt*x1k(7)) + (Jn*x1k(bp+k));
        end

    elseif sit==11 %GRE
        if sls(k)<33
            Ck(s,1) = rho + pc + (Jw*x1k(4)) + (Jt*x1k(5));
            Ck(f,1) = rho + lc + (Jw*x1k(4)) + (Jt*x1k(5)) + (Jn*x1k(bp+k));
        elseif sls(k)<59
            Ck(s,1) = rho + pc + (Jw*x1k(4)) + (Jt*x1k(6));
            Ck(f,1) = rho + lc + (Jw*x1k(4)) + (Jt*x1k(6)) + (Jn*x1k(bp+k));
        else
            Ck(s,1) = rho + pc + (Jw*x1k(4)) + (Jt*x1k(7));
            Ck(f,1) = rho + lc + (Jw*x1k(4)) + (Jt*x1k(7)) + (Jn*x1k(bp+k));
        end

    elseif sit==12 %GRC
        if sls(k)<33
            Ck(s,1) = rho + pc + (Jw*x1k(4)) + (Jt*x1k(5));
            Ck(f,1) = rho + lc + (Jw*x1k(4)) + (Jt*x1k(5)) + (Jn*x1k(bp+k));
        elseif sls(k)<59
            Ck(s,1) = rho + pc + (Jw*x1k(4)) + (Jt*x1k(6));
            Ck(f,1) = rho + lc + (Jw*x1k(4)) + (Jt*x1k(6)) + (Jn*x1k(bp+k));
        elseif sls(k)<111
            Ck(s,1) = rho + pc + (Jw*x1k(4)) + (Jt*x1k(7));
            Ck(f,1) = rho + lc + (Jw*x1k(4)) + (Jt*x1k(7)) + (Jn*x1k(bp+k));
        else
            Ck(s,1) = rho + pc + (Jw*x1k(4)) + (Jt*x1k(8));
            Ck(f,1) = rho + lc + (Jw*x1k(4)) + (Jt*x1k(8)) + (Jn*x1k(bp+k));
        end
     elseif sit==13 %GEC
        if sls(k)<33
            Ck(s,1) = rho + pc + (Jw*x1k(4)) + (Jt*x1k(5));
            Ck(f,1) = rho + lc + (Jw*x1k(4)) + (Jt*x1k(5)) + (Jn*x1k(bp+k));
        elseif sls(k)<95
            Ck(s,1) = rho + pc + (Jw*x1k(4)) + (Jt*x1k(6));
            Ck(f,1) = rho + lc + (Jw*x1k(4)) + (Jt*x1k(6)) + (Jn*x1k(bp+k));
        elseif sls(k)<111
            Ck(s,1) = rho + pc + (Jw*x1k(4)) + (Jt*x1k(7));
            Ck(f,1) = rho + lc + (Jw*x1k(4)) + (Jt*x1k(7)) + (Jn*x1k(bp+k));
        else
            Ck(s,1) = rho + pc + (Jw*x1k(4)) + (Jt*x1k(8));
            Ck(f,1) = rho + lc + (Jw*x1k(4)) + (Jt*x1k(8)) + (Jn*x1k(bp+k));
        end

    elseif sit==14 %REC
        if sls(k)<59
            Ck(s,1) = rho + pc + (Jw*x1k(4)) + (Jt*x1k(5));
            Ck(f,1) = rho + lc + (Jw*x1k(4)) + (Jt*x1k(5)) + (Jn*x1k(bp+k));
        elseif sls(k)<95
            Ck(s,1) = rho + pc + (Jw*x1k(4)) + (Jt*x1k(6));
            Ck(f,1) = rho + lc + (Jw*x1k(4)) + (Jt*x1k(6)) + (Jn*x1k(bp+k));
        elseif sls(k)<111
            Ck(s,1) = rho + pc + (Jw*x1k(4)) + (Jt*x1k(7));
            Ck(f,1) = rho + lc + (Jw*x1k(4)) + (Jt*x1k(7)) + (Jn*x1k(bp+k));
        else
            Ck(s,1) = rho + pc + (Jw*x1k(4)) + (Jt*x1k(8));
            Ck(f,1) = rho + lc + (Jw*x1k(4)) + (Jt*x1k(8)) + (Jn*x1k(bp+k));
        end
    elseif sit==15 %GREC
        if sls(k)<33
            Ck(s,1) = rho + pc + (Jw*x1k(4)) + (Jt*x1k(5));
            Ck(f,1) = rho + lc + (Jw*x1k(4)) + (Jt*x1k(5)) + (Jn*x1k(bp+k));
        elseif sls(k)<59
            Ck(s,1) = rho + pc + (Jw*x1k(4)) + (Jt*x1k(6));
            Ck(f,1) = rho + lc + (Jw*x1k(4)) + (Jt*x1k(6)) + (Jn*x1k(bp+k));
        elseif sls(k)<95
            Ck(s,1) = rho + pc + (Jw*x1k(4)) + (Jt*x1k(7));
            Ck(f,1) = rho + lc + (Jw*x1k(4)) + (Jt*x1k(7)) + (Jn*x1k(bp+k));
        elseif sls(k)<111
            Ck(s,1) = rho + pc + (Jw*x1k(4)) + (Jt*x1k(8));
            Ck(f,1) = rho + lc + (Jw*x1k(4)) + (Jt*x1k(8)) + (Jn*x1k(bp+k));
        else
            Ck(s,1) = rho + pc + (Jw*x1k(4)) + (Jt*x1k(9));
            Ck(f,1) = rho + lc + (Jw*x1k(4)) + (Jt*x1k(9)) + (Jn*x1k(bp+k));
        end
    end

    % for ambiguity resolution
    Rw(s,1) = rho + pc;
    Rw(f,1) = rho + lc;
    
    ab1 = (f1^2)/(f1^2 - f2^2);
    ab2 = (f2^2)/(f1^2 - f2^2);
    elv = meas((4*k - 1),26);
    if sls(k)<33
        kp = options.vrt.gps.p; 
        kl = options.vrt.gps.l;
        
    elseif sls(k)<59
        kp = options.vrt.glo.p; 
        kl = options.vrt.glo.l;

    elseif sls(k)<95
        kp = options.vrt.gal.p;
        kl = options.vrt.gal.l;

    elseif sls(k)<111
        lss = [95,96,97,98,99];
        if any(lss==sls(k))
            kp = 10; 
            kl = 10;
        else
            kp = options.vrt.bds.p;
            kl = options.vrt.bds.l;
        end
    else
        lss = [153,154,155];
        if any(lss==sls(k))
            kp = 10; 
            kl = 10;
        else
            kp = options.vrt.bds.p;
            kl = options.vrt.bds.l;
        end

    end

    code_var = ((options.CodeStd *kp)^2)*(ab1+ab2);
    phas_var = ((options.PhaseStd*kl)^2)*(ab1+ab2);
    
    Rk(s,s) = code_var/(sind(elv)^2);
    Rk(f,f) = phas_var/(sind(elv)^2);
end

if any(isnan(Zk)) || any(isnan(Ck))
    [a,~] = find(isnan(Zk));
    for b=a
        Zk(b,:) = [];
        Ck(b,:) = [];
        Hk(b,:) = [];
        Rk(b,:) = [];
        Rk(:,b) = [];
    end
end

sres0 = zeros(size(Zk,1),1);
nt = 0;
while 1
    nt = nt + 1;
    Vk   = Zk - Ck;
    Sk = (Hk*p1k*Hk') + Rk;
    Kk = ((p1k*Hk'))/(Sk);

    DD = pinv(Hk'*Hk);
    kof = DD(1:bp,1:bp);

    dx   = Kk*Vk;                 
    xk = x1k + dx;                
    tnk = (eye(pno) - Kk*Hk);
    pk = tnk*p1k*tnk' + Kk*Rk*Kk';

    res  = (Ck + (Hk*dx)) - Zk;   
    vres = Rk - (Hk*pk*Hk');
    
    s0 = (Vk'*(Sk\Vk))/(size(Vk,1));

    sres = abs(res)./sqrt(s0.*abs(diag(vres)));
    dres = abs(sres - sres0);
    if any(dres>0.5)
        mm = find(sres == max(sres));
        k0 = 2.5; k1 = 8;
        if sres(mm,1)>k1
            sm = 1*10^-10;
            Rk(mm,mm) = Rk(mm,mm)/sm;
        elseif sres(mm,1)>k0
            hp = abs(sres(mm,1));
            sm = (k0/hp)*((k1 - hp)/(k1 - k0))^2;
            Rk(mm,mm) = Rk(mm,mm)/sm;
        end
        sres0 = sres;
    else    
        res  = (Ck + (Hk*dx)) - Zk;
        break
    end
end

% 
amb.hk = Hk;
amb.rk = Rk;
amb.vk = Zk - (Rw + (Hk(:,4:bp)*xk(4:bp)));
amb.xk = xk;
amb.pk = pk;
amb.x1 = x1k;

end