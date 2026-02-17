function [xs,pks,kofs,ress,satno,fxpos,tfxl] = modelandfiltering(data,options,waitbarHandle,progressStart,progressEnd)

% Check if waitbar parameters are provided
useWaitbar = false;
if nargin >= 3 && ~isempty(waitbarHandle)
    useWaitbar = true;
    if nargin < 4
        progressStart = 0.5; % Default start progress
    end
    if nargin < 5
        progressEnd = 0.9; % Default end progress
    end
end

rn = size(data.rnx.pob,1);
cn = size(data.rnx.pob,2);

% Update waitbar with initial message if available
if useWaitbar
    waitbar(progressStart, waitbarHandle, {'Please wait...' 'Initializing modeling process...'});
end

cv = 299792458;
[freq,wavl] = frequencies;

mopt = zeros(1,11);
mopt( 1) = options.SatClk;
mopt( 2) = options.SatAPC;
mopt( 3) = options.RecAPC;
mopt( 4) = options.RecARP;
mopt( 5) = options.RelClk;
mopt( 6) = options.SatWind;
mopt( 7) = options.AtmTrop;
mopt( 9) = options.RelPath;
mopt(11) = options.Solid;

if strcmp(data.rnx.head.time.system,'GPS') || isempty(data.rnx.head.time.system)
    dt = 51.184;
else
    frst = data.rnx.head.time.first;
    [~,mjd] = cal2jul(frst(1), frst(2), frst(3), (frst(4)*3600 + frst(5)*60 + frst(6)));
    leap = leapsec(mjd);
    dt  = 32.184 + leap;
end

arp = data.rnx.head.ant.hen;
if size(arp,1)~=1
    arp = arp';
end

atx = data.atx;
rapc = data.atx.rcv;

[bp,~,sit] = d_system(data);

fs = bp + cn;
xs   = zeros(rn,fs);
pks  = zeros(fs,fs,rn);
kofs = zeros(bp,bp,rn);
ress = NaN(rn,cn,1);
satno = NaN(rn,1);
pres = NaN(rn,cn,2);

refs = NaN(1,5);
stes = NaN(1,5);
sdwl = NaN(rn,cn);
tfxl = NaN(rn,cn,2);

fxll  = NaN(1,cn);
fxpos = NaN(rn,3);

% Calculate progress increment per iteration if waitbar is being used
if useWaitbar
    progressRange = progressEnd - progressStart;
    progressStep = progressRange / rn;
end

for i=1:rn
    % Update waitbar with current epoch progress
    if useWaitbar
        currentProgress = progressStart + ((i-1) * progressStep);
        waitbar(currentProgress, waitbarHandle, {'Modeling and filtering...' sprintf('Modeling epoch %d of %d (%.2f%%)', i, rn, currentProgress*100)});
    end
    
    year = data.rnx.head.time.first(1); mon = data.rnx.head.time.first(2); day = data.rnx.head.time.first(3);
    doy = data.rnx.head.time.doy; 
    secod = data.rnx.obs.ep(i,1);
    [~,mjd_tt] = cal2jul(year,mon,day,(secod+dt));        
    [~,mjd]    = cal2jul(year,mon,day,(secod));

    sun_xyz = c_sunpos(mjd_tt); 
    mon_xyz = c_moonpos(mjd_tt);

    r_xyz = data.rnx.head.rec.pos;
    if size(r_xyz,1)~=1
        r_xyz = r_xyz';
    end

    [ellp] = xyz2plh(r_xyz,0);
    dlat = ellp(1);
    dlon = ellp(2);
    hell = ellp(3);

    grid = load('gpt3_5.mat');
    [p,~,~,Tm,e,ah,aw,la,~,~,~,~,~] = gpt3_5_fast(mjd,dlat,dlon,hell,0,grid.gpt3_grid);

    sats = find(data.rnx.pob(i,:) == 1);
    mn = size(sats,2)*4;
    model  = zeros(mn,32);

    t = 0;
    for k = sats
        if k<33
            sno = k;
            sapc = atx.gps;
            kc = 1;
        elseif k<59
            sno = k-32;
            sapc = atx.glo;
            kc = 2;
        elseif k<95
            sno = k-58;
            sapc = atx.gal;
            kc = 3;
        else
            sno = k-94;
            sapc = atx.bds;
            kc = 4;
        end

        for u=1:4
            t = t + 1;
            
            model(t,1) = year;
            model(t,2) = doy;
            model(t,3) = secod;
            model(t,4) = k;
            model(t,5) = data.sat.tfs(i,k);

            s_xyz = data.sat.crt(i,1:3,k);
            v_xyz = data.sat.crt(i,4:6,k);
            azi   = data.sat.lcl(i,k,1);
            elv   = data.sat.lcl(i,k,2);
            
            switch u
                case 1 
                    model(t,6 ) = data.rnx.sob(i,k,1);

                    sdt = c_satapc(s_xyz,r_xyz,sun_xyz,sapc,sno,1);
                    if ~isnan(sdt), model(t,16) = sdt; end

                    kn = 2*kc - 1;
                    model(t,17) = c_recapc(s_xyz,r_xyz,rapc,kn,azi,elv);

                case 2
                    model(t,6) = data.rnx.sob(i,k,3);
                    
                    sdt = c_satapc(s_xyz,r_xyz,sun_xyz,sapc,sno,2);
                    if ~isnan(sdt), model(t,16) = sdt; end

                    kn = 2*kc;
                    model(t,17) = c_recapc(s_xyz,r_xyz,rapc,kn,azi,elv);
                case 3
                    model(t,6) = data.rnx.sob(i,k,2);
                    
                    sdt = c_satapc(s_xyz,r_xyz,sun_xyz,sapc,sno,1);
                    if ~isnan(sdt), model(t,16) = sdt; end
                    
                    kn = 2*kc - 1;
                    model(t,17) = c_recapc(s_xyz,r_xyz,rapc,kn,azi,elv);
                    
                    if i==1
                        prev=0;
                    else
                        if ~isnan(pres(i-1,k,1))
                            prev = pres(i-1,k,1);
                        else
                            prev = 0;
                        end
                    end
                    wu = c_windup(r_xyz,s_xyz,sun_xyz,prev);
                    pres(i,k,1) = wu;
                    model(t,20) = wu*wavl(k,1);
                    
                case 4
                    model(t,6) = data.rnx.sob(i,k,4);
                    sdt = c_satapc(s_xyz,r_xyz,sun_xyz,sapc,sno,2);
                    if ~isnan(sdt), model(t,16) = sdt; end
                    
                    kn = 2*kc;
                    model(t,17) = c_recapc(s_xyz,r_xyz,rapc,kn,azi,elv);
                    
                    if i==1
                        prev=0;
                    else
                        if ~isnan(pres(i-1,k,2))
                            prev = pres(i-1,k,2);
                        else
                            prev = 0;
                        end
                    end
                    wu = c_windup(r_xyz,s_xyz,sun_xyz,prev);
                    pres(i,k,2) = wu;
                    model(t,20) = wu*wavl(k,2); 
            end
        
            model(t,8:10)  = s_xyz;
            
            model(t,11:13) = v_xyz;
            
            model(t,14) = norm(s_xyz - r_xyz);
            
            model(t,15) = -(data.sat.crt(i,7,k)*cv);
            
            model(t,18) = c_recarp(s_xyz,r_xyz,arp);
            
            model(t,19) = -(c_relclk(s_xyz,v_xyz));

            [Trop,Mwet,Mn,Me,ZHD,ZWD] = c_trop(r_xyz,s_xyz,mjd,p,ah,aw,Tm,e,la);
            
            model(t,21) = Trop;
            model(t,28) = Mwet;
            model(t,29) = Mn;
            model(t,30) = Me;
            model(t,31) = ZHD;
            model(t,32) = ZWD;
            
            model(t,23) = c_rpath(r_xyz,s_xyz);
            
            model(t,25) = c_solid(r_xyz,s_xyz,sun_xyz,mon_xyz);
            
            model(t,26) = elv;
            
            model(t,27) = azi;
            
            full = model(t,15:25);
            model(t,7) = sum(full(mopt==1));
        end
    end
    
    sls = find(data.rnx.pob(i,:) == 1);
    sno = size(sls,2);
    pno = bp + sno;
    meas = model;
    Q  = zeros(pno);

    tint = data.rnx.head.time.int;
    
    if options.ProMod == 1
        Q(1,1) = (options.NosPos*10^(options.NosPos2))*(tint);
        Q(2,2) = (options.NosPos*10^(options.NosPos2))*(tint);
        Q(3,3) = (options.NosPos*10^(options.NosPos2))*(tint);
    end
    Q(4,4) = (options.NosTrop*10^(options.NosTrop2))*(tint);
    Q(5,5) = (options.NosClk*10^(options.NosClk2))*(tint);

    if bp==6
        Q(6,6) = (options.NosClk*10^(options.NosClk2))*(tint);
    elseif bp==7
        Q(6,6) = (options.NosClk*10^(options.NosClk2))*(tint);
        Q(7,7) = (options.NosClk*10^(options.NosClk2))*(tint);
    elseif bp==8
        Q(6,6) = (options.NosClk*10^(options.NosClk2))*(tint);
        Q(7,7) = (options.NosClk*10^(options.NosClk2))*(tint);
        Q(8,8) = (options.NosClk*10^(options.NosClk2))*(tint);
    elseif bp==9
        Q(6,6) = (options.NosClk*10^(options.NosClk2))*(tint);
        Q(7,7) = (options.NosClk*10^(options.NosClk2))*(tint);
        Q(8,8) = (options.NosClk*10^(options.NosClk2))*(tint);
        Q(9,9) = (options.NosClk*10^(options.NosClk2))*(tint);
    end
    
    F = eye(pno);

    if i == 1
        x1k = zeros(pno,1);
        x1k(1:3,1) = data.rnx.head.rec.pos';

        p1k = zeros(pno);
        p1k(1,1) = (options.IntPos*10^(options.IntPos2))^2;
        p1k(2,2) = (options.IntPos*10^(options.IntPos2))^2;
        p1k(3,3) = (options.IntPos*10^(options.IntPos2))^2;
        p1k(4,4) = (options.IntTrop*10^(options.IntTrop2))^2;
        p1k(5,5) = (options.IntClk*10^(options.IntClk2))^2;
        if bp==6
            p1k(6,6) = (options.IntClk*10^(options.IntClk2))^2;
        elseif bp==7
            p1k(6,6) = (options.IntClk*10^(options.IntClk2))^2;
            p1k(7,7) = (options.IntClk*10^(options.IntClk2))^2;
        elseif bp==8
            p1k(6,6) = (options.IntClk*10^(options.IntClk2))^2;
            p1k(7,7) = (options.IntClk*10^(options.IntClk2))^2;
            p1k(8,8) = (options.IntClk*10^(options.IntClk2))^2;
        elseif bp==9
            p1k(6,6) = (options.IntClk*10^(options.IntClk2))^2;
            p1k(7,7) = (options.IntClk*10^(options.IntClk2))^2;
            p1k(8,8) = (options.IntClk*10^(options.IntClk2))^2;
            p1k(9,9) = (options.IntClk*10^(options.IntClk2))^2;
        end

        for u=(bp+1):pno
            p1k(u,u) = (options.IntAmb*10^(options.IntAmb2))^2;
        end
        
    else
        x1k = zeros(pno,1);
        x1k(1:bp,1) = xs(i-1,1:bp);
        for k=1:sno
            snm = sls(k);
            x1k(bp+k,1) = xs(i-1,bp+snm);
        end
        x1k = F*x1k;
        
        p1k = zeros(pno);
        for r=1:size(p1k,1)
            for c=1:size(p1k,2)
                if r<(bp+1) && c<(bp+1)
                    p1k(r,c) = ps(r,c);
                elseif r<(bp+1) && c>bp
                    sn = sls(c-bp);
                    p1k(r,c) = ps(r,sn+bp);
                elseif r>bp && c<(bp+1)
                    sn = sls(r-bp);
                    p1k(r,c) = ps(sn+bp,c);
                else
                    f1 = sls(r-bp);
                    f2 = sls(c-bp);
                    p1k(r,c) = ps(f1+bp,f2+bp);
                end
            end
        end
        
        for k=1:sno
            snm = bp + k;
            if p1k(snm,snm)==0
                p1k(snm,snm) = (options.IntAmb*10^(options.IntAmb2))^2;
            end
        end
        
        p1k = F*p1k*F' + Q;

        %[xk,pk,kof,res,amb] = f_rkalman(sno,sls,pno,meas,x1k,p1k,sit,bp,freq,options);
    end

    % Update waitbar with Kalman filter status
    if useWaitbar && mod(i, max(1, floor(rn/20))) == 0
        currentProgress = progressStart + ((i-0.5) * progressStep);
        waitbar(currentProgress, waitbarHandle, {'Please wait...' sprintf('Kalman filtering epoch %d of %d', i, rn)});
    end

    [xk,pk,kof,res,amb] = f_rkalman(sno,sls,pno,meas,x1k,p1k,sit,bp,freq,options);
    
    [refs,stes,sdwl,fxwl,fxll,fpos] = f_ambiguity3(refs,meas,data,stes,i,sdwl,bp,sls,amb,fxll);

    tfxl(i,:,1) = fxwl;
    tfxl(i,:,2) = fxll;
    fxpos(i,:) = fpos;

    xs(i,1:bp) = xk(1:bp,1);
    for k = 1:sno
       snm = sls(k) + bp;
       xs(i,snm)  = xk(bp+k,1);
    end

    ps = zeros(fs);
    for r=1:size(pk,1)
       for c=1:size(pk,2)
           if r<(bp+1) && c<(bp+1)
               ps(r,c) = pk(r,c);
           elseif r<(bp+1) && c>bp
               sn = sls(c-bp);
               ps(r,sn+bp) = pk(r,c);
           elseif r>bp && c<(bp+1)
               sn = sls(r-bp);
               ps(sn+bp,c) = pk(r,c);
           else
               sn1 = sls(r-bp); sn2 = sls(c-bp);
               ps(sn1+bp,sn2+bp) = pk(r,c);
           end
       end
    end
    pks(:,:,i) = ps;
    
    kofs(:,:,i) = kof(:,:,1);
    
    for rr=1:round(size(res,1)/2)
        sns = sls(rr);
        sts = (rr-1)*2 + 1;
        fns = 2*rr;
        ress(i,sns,1) = res(sts);
        ress(i,sns,2) = res(fns);
    end
    
    sns = find(data.rnx.pob(i,:) == 1);
    satno(i) = size(sns,2); 
    
    if options.ProMod == 1
        data.rnx.head.rec.pos = xk(1:3)';
    end
    
    % Update waitbar with completion percentage
    if useWaitbar
        currentProgress = progressStart + (i * progressStep);
        waitbar(currentProgress, waitbarHandle, {'Modeling and filtering...' sprintf('Completed epoch %d of %d (%.2f%%)', i, rn, currentProgress*100)});
    end
end

% Final waitbar update
if useWaitbar
    waitbar(progressEnd, waitbarHandle, {'Please wait...' 'Modeling and filtering completed'});
end

end