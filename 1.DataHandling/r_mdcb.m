function [dcb] = r_mdcb(fname)

dcb.gps.c1cc1w = zeros(32,1);
dcb.gps.c1wc2w = zeros(32,1);
dcb.glo.c1cc1p = zeros(26,1);
dcb.glo.c1pc2p = zeros(26,1);
dcb.gal.c1cc5q = zeros(36,1);
dcb.bds.c2ic6i = zeros(61,1);
dcb.bds.c2ic7i = zeros(61,1);

fid = fopen(fname);
while ~feof(fid)
    line = fgetl(fid);
    % GPS
    if strcmp(strtrim(line(1:4)),'DSB') && strcmp(line(12),'G')
        k = sscanf(line(13:14),'%d');
        if strcmp(line(26:33),'C1C  C1W')
            tm = sscanf(line(85:91),'%f');
            dcb.gps.c1cc1w(k,1) = tm*10^-9;%seconds
        end
        
        if strcmp(line(26:33),'C1W  C2W')
            tm = sscanf(line(85:91),'%f');
            dcb.gps.c1wc2w(k,1) = tm*10^-9;%seconds
        end
    % GLONASS
    elseif strcmp(strtrim(line(1:4)),'DSB') && strcmp(line(12),'R')
        k = sscanf(line(13:14),'%d');
        if strcmp(line(26:33),'C1C  C1P')
            tm = sscanf(line(85:91),'%f');
            dcb.glo.c1cc1p(k,1) = tm*10^-9;
        end
        
        if strcmp(line(26:33),'C1P  C2P')
            tm = sscanf(line(85:91),'%f');
            dcb.glo.c1pc2p(k,1) = tm*10^-9;
        end
    % Galileo
    elseif strcmp(strtrim(line(1:4)),'DSB') && strcmp(line(12),'E')
        k = sscanf(line(13:14),'%d');
        if strcmp(line(26:33),'C1C  C5Q')
            tm = sscanf(line(85:91),'%f');
            dcb.gal.c1cc5q(k,1) = tm*10^-9;
        end
    % BDS
    elseif strcmp(strtrim(line(1:4)),'DSB') && strcmp(line(12),'C')
        k = sscanf(line(13:14),'%d');
        if strcmp(line(26:33),'C2I  C6I')
            tm = sscanf(line(85:91),'%f');
            dcb.bds.c2ic6i(k,1) = tm*10^-9;
        end
        if strcmp(line(26:33),'C2I  C7I')
            tm = sscanf(line(85:91),'%f');
            dcb.bds.c2ic7i(k,1) = tm*10^-9;
        end
    end
end

fclose(fid);
end