function [data] = r_osb(fname,data)

% GPS code
osb.gps.c1c = zeros(32,1); osb.gps.c1l = zeros(32,1);
osb.gps.c1x = zeros(32,1); osb.gps.c1w = zeros(32,1);
osb.gps.c2l = zeros(32,1); osb.gps.c2s = zeros(32,1);
osb.gps.c2x = zeros(32,1); osb.gps.c2w = zeros(32,1);
% GPS phase
osb.gps.l1c = zeros(32,1); osb.gps.l1l = zeros(32,1);
osb.gps.l1x = zeros(32,1); osb.gps.l1w = zeros(32,1);
osb.gps.l2l = zeros(32,1); osb.gps.l2s = zeros(32,1);
osb.gps.l2x = zeros(32,1); osb.gps.l2w = zeros(32,1);

% GLONASS code
osb.glo.c1p = zeros(26,1); osb.glo.c2p = zeros(26,1);
% GLONASS phase
osb.glo.l1p = zeros(26,1); osb.glo.l2p = zeros(26,1);

% Galileo code
osb.gal.c1c = zeros(36,1); osb.gal.c1x = zeros(36,1);
osb.gal.c5i = zeros(36,1); osb.gal.c5q = zeros(36,1);
osb.gal.c5x = zeros(36,1);
% Galileo phase
osb.gal.l1c = zeros(36,1); osb.gal.l1x = zeros(36,1); 
osb.gal.l5i = zeros(36,1); osb.gal.l5q = zeros(36,1);
osb.gal.l5x = zeros(36,1);
% BeiDou code
osb.bds.c2i = zeros(62,1); osb.bds.c6i = zeros(62,1);
% BeiDou phase
osb.bds.l2i = zeros(62,1); osb.bds.l6i = zeros(62,1);

try
    fid = fopen(fname);

    while ~feof(fid)
        line = fgetl(fid);
        if length(line)>12
            if strcmp(strtrim(line(1:4)),'OSB') && strcmp(line(12),'G') ...
                    && isempty(strip(line(16:24)))
                sno = sscanf(line(13:14),'%d');
                sgn = line(26:28);
                tm = sscanf(line(71:91),'%f')*10^-9;
                if strcmp(sgn,'C1C')
                    osb.gps.c1c(sno,1) = tm;
                elseif strcmp(sgn,'C1L')
                    osb.gps.c1l(sno,1) = tm;
                elseif strcmp(sgn,'C1X')
                    osb.gps.c1x(sno,1) = tm;
                elseif strcmp(sgn,'C1W')
                    osb.gps.c1w(sno,1) = tm;
                elseif strcmp(sgn,'C2L')
                    osb.gps.c2l(sno,1) = tm;
                elseif strcmp(sgn,'C2S')
                    osb.gps.c2s(sno,1) = tm;
                elseif strcmp(sgn,'C2X')
                    osb.gps.c2x(sno,1) = tm;
                elseif strcmp(sgn,'C2W')
                    osb.gps.c2w(sno,1) = tm;
                elseif strcmp(sgn,'L1C')
                    osb.gps.l1c(sno,1) = tm;
                elseif strcmp(sgn,'L1L')
                    osb.gps.l1l(sno,1) = tm;
                elseif strcmp(sgn,'L1X')
                    osb.gps.l1x(sno,1) = tm;
                elseif strcmp(sgn,'L1W')
                    osb.gps.l1w(sno,1) = tm;
                elseif strcmp(sgn,'L2L')
                    osb.gps.l2l(sno,1) = tm;
                elseif strcmp(sgn,'L2S')
                    osb.gps.l2s(sno,1) = tm;
                elseif strcmp(sgn,'L2X')
                    osb.gps.l2x(sno,1) = tm;
                elseif strcmp(sgn,'L2W')
                    osb.gps.l2w(sno,1) = tm; 
                end
    
            elseif strcmp(strtrim(line(1:4)),'OSB') && strcmp(line(12),'R')...
                    && isempty(strip(line(16:24)))
                sno = sscanf(line(13:14),'%d');
                sgn = line(26:28);
                tm = sscanf(line(71:91),'%f')*10^-9;
                if strcmp(sgn,'C1P')
                    osb.glo.c1p(sno,1) = tm;
                elseif strcmp(sgn,'C2P')
                    osb.glo.c2p(sno,1) = tm;
                elseif strcmp(sgn,'L1P')
                    osb.glo.l1p(sno,1) = tm;
                elseif strcmp(sgn,'L2P')
                    osb.glo.l2p(sno,1) = tm;
                end
    
            elseif strcmp(strtrim(line(1:4)),'OSB') && strcmp(line(12),'E')...
                    && isempty(strip(line(16:24)))
                sno = sscanf(line(13:14),'%d');
                sgn = line(26:28);
                tm = sscanf(line(71:91),'%f')*10^-9;
                if strcmp(sgn,'C1C')
                    osb.gal.c1c(sno,1) = tm;
                elseif strcmp(sgn,'C1X')
                    osb.gal.c1x(sno,1) = tm;
                elseif strcmp(sgn,'C5I')
                    osb.gal.c5i(sno,1) = tm;
                elseif strcmp(sgn,'C5Q')
                    osb.gal.c5q(sno,1) = tm;
                elseif strcmp(sgn,'C5X')
                    osb.gal.c5x(sno,1) = tm;
                elseif strcmp(sgn,'L1C')
                    osb.gal.l1c(sno,1) = tm;
                elseif strcmp(sgn,'L1X')
                    osb.gal.l1x(sno,1) = tm;
                elseif strcmp(sgn,'L5I')
                    osb.gal.l5i(sno,1) = tm;
                elseif strcmp(sgn,'L5Q')
                    osb.gal.l5q(sno,1) = tm;
                elseif strcmp(sgn,'L5X')
                    osb.gal.l5x(sno,1) = tm;
                end
    
            elseif strcmp(strtrim(line(1:4)),'OSB') && strcmp(line(12),'C')...
                    && isempty(strip(line(16:24)))
                sno = sscanf(line(13:14),'%d');
                sgn = line(26:28);
                tm = sscanf(line(71:91),'%f')*10^-9;
                if strcmp(sgn,'C2I')
                    osb.bds.c2i(sno,1) = tm;
                elseif strcmp(sgn,'C6I')
                    osb.bds.c6i(sno,1) = tm;
                elseif strcmp(sgn,'L2I')
                    osb.bds.l2i(sno,1) = tm;
                elseif strcmp(sgn,'L6I')
                    osb.bds.l6i(sno,1) = tm;
                end
            end
        end
    end
    
    data.osb = osb;
    fclose(fid);
catch ME
    data.err.cod = 1;
    msg{1} = ['Error message:',ME.message];
    data.err.msg = msg;
    return
end

end