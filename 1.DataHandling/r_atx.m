function [data] = r_atx(fname,data)

if data.err.cod == 1
    return
end

try
    [fid,~] = fopen(fname);

    type = data.rnx.head.ant.type;

    atx.gps.daz = zeros(32,4);
    atx.glo.daz = zeros(26,4);
    atx.gal.daz = zeros(36,4);
    atx.bds.daz = zeros(62,4);

    atx.gps.neu = zeros(32,3,2);
    atx.glo.neu = zeros(26,3,2);
    atx.gal.neu = zeros(36,3,2);
    atx.bds.neu = zeros(62,3,2);

    atx.gps.nad = cell(32,2);
    atx.glo.nad = cell(26,2);
    atx.gal.nad = cell(36,2);
    atx.bds.nad = cell(62,2);

    atx.gps.adp = cell(32,2);
    atx.glo.adp = cell(26,2);
    atx.gal.adp = cell(36,2);
    atx.bds.adp = cell(62,2);

    atx.rcv.neu = zeros(1,3,8);
    atx.rcv.nad = cell(1,8);
    atx.rcv.adp = cell(1,8);
    atx.rcv.daz = zeros(1,4);
    
    linenum = 0;
    while ~feof(fid)
        tline = fgetl(fid); linenum = linenum + 1;
        tag   = strtrim(tline(61:end));
        if strcmp(tag,'START OF ANTENNA')
            tline = fgetl(fid); linenum = linenum + 1;
            tag   = strtrim(tline(61:end));

            if strcmp(tag,'TYPE / SERIAL NO') && strcmp(tline(21),'G')
                sno = sscanf(tline(22:23),'%d');
                if sno>size(atx.gps.neu,1)
                    continue
                end
                
                while ~strcmp(tag,'END OF ANTENNA')
                    tline = fgetl(fid); linenum = linenum + 1;
                    tag   = strtrim(tline(61:end));
    
                    if strcmp(tag,'DAZI')
                        dazi = sscanf(tline(3:8),'%f');
                        atx.gps.daz(sno,1) = dazi;
                    end
                    if strcmp(tag,'ZEN1 / ZEN2 / DZEN')
                        dzen = sscanf(tline(3:20),'%f');
                        atx.gps.daz(sno,2:4) = dzen';
                        nzen = (dzen(2)-dzen(1))/dzen(3) + 1;
                    end
                    if strcmp(tag,'START OF FREQUENCY') && strcmp(tline(4:6),'G01')
                        fno = 1;
                        tline = fgetl(fid); linenum = linenum + 1;
                        tag   = strtrim(tline(61:end));
                        if strcmp(tag,'NORTH / EAST / UP')
                            atx.gps.neu(sno,:,fno) = sscanf(tline,'%f',[1,3]);
                            tline = fgetl(fid); linenum = linenum + 1;
                            atx.gps.nad{sno,fno} = sscanf(tline(9:end),'%f',[1,nzen]);
                            if dazi~=0
                                nn = 360/dazi + 1;
                                adp = NaN(nn,nzen+1);
                                for i=1:nn
                                    tline = fgetl(fid);
                                    linenum = linenum + 1;
                                    adp(i,:) = sscanf(tline,'%f',[1,nzen+1]);
                                end
                                atx.gps.adp{sno,fno} = adp;
                            end
                        end
                    elseif strcmp(tag,'START OF FREQUENCY') && strcmp(tline(4:6),'G02')
                        fno = 2;
                        tline = fgetl(fid); linenum = linenum + 1;
                        tag   = strtrim(tline(61:end));
                        if strcmp(tag,'NORTH / EAST / UP')
                            atx.gps.neu(sno,:,fno) = sscanf(tline,'%f',[1,3]);
                            tline = fgetl(fid); linenum = linenum + 1;
                            atx.gps.nad{sno,fno} = sscanf(tline(9:end),'%f',[1,nzen]);
                            if dazi~=0
                                nn = 360/dazi + 1;
                                adp = NaN(nn,nzen+1);
                                for i=1:nn
                                    tline = fgetl(fid);
                                    linenum = linenum + 1;
                                    adp(i,:) = sscanf(tline,'%f',[1,nzen+1]);
                                end
                                atx.gps.adp{sno,fno} = adp;
                            end
                        end
                    end
                end

            elseif strcmp(tag,'TYPE / SERIAL NO') && strcmp(tline(21),'R')
                sno = sscanf(tline(22:23),'%d');
                if sno>size(atx.glo.neu,1) 
                    continue
                end

                while ~strcmp(tag,'END OF ANTENNA')
                    tline = fgetl(fid); linenum = linenum + 1;
                    tag   = strtrim(tline(61:end));
    
                    if strcmp(tag,'DAZI')
                        dazi = sscanf(tline(3:8),'%f');
                        atx.glo.daz(sno,1) = dazi;
                    end
                    if strcmp(tag,'ZEN1 / ZEN2 / DZEN')
                        dzen = sscanf(tline(3:20),'%f');
                        atx.glo.daz(sno,2:4) = dzen';
                        nzen = (dzen(2)-dzen(1))/dzen(3) + 1;
                    end
                    if strcmp(tag,'START OF FREQUENCY') && strcmp(tline(4:6),'R01')
                        fno = 1;
                        tline = fgetl(fid); linenum = linenum + 1;
                        tag   = strtrim(tline(61:end));
                        if strcmp(tag,'NORTH / EAST / UP')
                            atx.glo.neu(sno,:,fno) = sscanf(tline,'%f',[1,3]);
                            tline = fgetl(fid); linenum = linenum + 1;
                            atx.glo.nad{sno,fno} = sscanf(tline(9:end),'%f',[1,nzen]);
                            if dazi~=0
                                nn = 360/dazi + 1;
                                adp = NaN(nn,nzen+1);
                                for i=1:nn
                                    tline = fgetl(fid);
                                    linenum = linenum + 1;
                                    adp(i,:) = sscanf(tline,'%f',[1,nzen+1]);
                                end
                                atx.glo.adp{sno,fno} = adp;
                            end
                        end
                    elseif strcmp(tag,'START OF FREQUENCY') && strcmp(tline(4:6),'R02')
                        fno = 2;
                        tline = fgetl(fid); linenum = linenum + 1;
                        tag   = strtrim(tline(61:end));
                        if strcmp(tag,'NORTH / EAST / UP')
                            atx.glo.neu(sno,:,fno) = sscanf(tline,'%f',[1,3]);
                            tline = fgetl(fid); linenum = linenum + 1;
                            atx.glo.nad{sno,fno} = sscanf(tline(9:end),'%f',[1,nzen]);
                            if dazi~=0
                                nn = 360/dazi + 1;
                                adp = NaN(nn,nzen+1);
                                for i=1:nn
                                    tline = fgetl(fid);
                                    linenum = linenum + 1;
                                    adp(i,:) = sscanf(tline,'%f',[1,nzen+1]);
                                end
                                atx.glo.adp{sno,fno} = adp;
                            end
                        end
                    end
                end

            elseif strcmp(tag,'TYPE / SERIAL NO') && strcmp(tline(21),'E')
                sno = sscanf(tline(22:23),'%d');
                if sno>size(atx.gal.neu,1)
                    continue
                end

                while ~strcmp(tag,'END OF ANTENNA')
                    tline = fgetl(fid); linenum = linenum + 1;
                    tag   = strtrim(tline(61:end));
    
                    if strcmp(tag,'DAZI')
                        dazi = sscanf(tline(3:8),'%f');
                        atx.gal.daz(sno,1) = dazi;
                    end 
                    if strcmp(tag,'ZEN1 / ZEN2 / DZEN')
                        dzen = sscanf(tline(3:20),'%f');
                        atx.gal.daz(sno,2:4) = dzen';
                        nzen = (dzen(2)-dzen(1))/dzen(3) + 1;
                    end
                    if strcmp(tag,'START OF FREQUENCY') && strcmp(tline(4:6),'E01')
                        fno = 1;
                        tline = fgetl(fid); linenum = linenum + 1;
                        tag   = strtrim(tline(61:end));
                        if strcmp(tag,'NORTH / EAST / UP')
                            atx.gal.neu(sno,:,fno) = sscanf(tline,'%f',[1,3]);
                            tline = fgetl(fid); linenum = linenum + 1;
                            atx.gal.nad{sno,fno} = sscanf(tline(9:end),'%f',[1,nzen]);
                            if dazi~=0
                                nn = 360/dazi + 1;
                                adp = NaN(nn,nzen+1);
                                for i=1:nn
                                    tline = fgetl(fid);
                                    linenum = linenum + 1;
                                    adp(i,:) = sscanf(tline,'%f',[1,nzen+1]);
                                end
                                atx.gal.adp{sno,fno} = adp;
                            end
                        end
                    elseif strcmp(tag,'START OF FREQUENCY') && strcmp(tline(4:6),'E05')
                        fno = 2;
                        tline = fgetl(fid); linenum = linenum + 1;
                        tag   = strtrim(tline(61:end));
                        if strcmp(tag,'NORTH / EAST / UP')
                            atx.gal.neu(sno,:,fno) = sscanf(tline,'%f',[1,3]);
                            tline = fgetl(fid); linenum = linenum + 1;
                            atx.gal.nad{sno,fno} = sscanf(tline(9:end),'%f',[1,nzen]);
                            if dazi~=0
                                nn = 360/dazi + 1;
                                adp = NaN(nn,nzen+1);
                                for i=1:nn
                                    tline = fgetl(fid);
                                    linenum = linenum + 1;
                                    adp(i,:) = sscanf(tline,'%f',[1,nzen+1]);
                                end
                                atx.gal.adp{sno,fno} = adp;
                            end
                        end
                    end
                end

            elseif strcmp(tag,'TYPE / SERIAL NO') && strcmp(tline(21),'C')
                sno = sscanf(tline(22:23),'%d');
                if sno>size(atx.bds.neu,1)
                    continue
                end
                while ~strcmp(tag,'END OF ANTENNA')
                    tline = fgetl(fid); linenum = linenum + 1;
                    tag   = strtrim(tline(61:end));
    
                    if strcmp(tag,'DAZI')
                        dazi = sscanf(tline(3:8),'%f');
                        atx.bds.daz(sno,1) = dazi;
                    end
                    if strcmp(tag,'ZEN1 / ZEN2 / DZEN')
                        dzen = sscanf(tline(3:20),'%f');
                        atx.bds.daz(sno,2:4) = dzen';
                        nzen = (dzen(2)-dzen(1))/dzen(3) + 1;
                    end
                    if strcmp(tag,'START OF FREQUENCY') && strcmp(tline(4:6),'C02')
                        fno = 1;
                        tline = fgetl(fid); linenum = linenum + 1;
                        tag   = strtrim(tline(61:end));
                        if strcmp(tag,'NORTH / EAST / UP')
                            atx.bds.neu(sno,:,fno) = sscanf(tline,'%f',[1,3]);
                            tline = fgetl(fid);
                            linenum = linenum + 1;
                            atx.bds.nad{sno,fno} = sscanf(tline(9:end),'%f',[1,nzen]);
                            if dazi~=0
                                nn = 360/dazi + 1;
                                adp = NaN(nn,nzen+1);
                                for i=1:nn
                                    tline = fgetl(fid);
                                    linenum = linenum + 1;
                                    adp(i,:) = sscanf(tline,'%f',[1,nzen+1]);
                                end
                                atx.bds.adp{sno,fno} = adp;
                            end
                        end
                    elseif strcmp(tag,'START OF FREQUENCY') && strcmp(tline(4:6),'C06')
                        fno = 2;
                        tline = fgetl(fid); linenum = linenum + 1;
                        tag   = strtrim(tline(61:end));
                        if strcmp(tag,'NORTH / EAST / UP')
                            atx.bds.neu(sno,:,fno) = sscanf(tline,'%f',[1,3]);
                            tline = fgetl(fid); linenum = linenum + 1;
                            atx.bds.nad{sno,fno} = sscanf(tline(9:end),'%f',[1,nzen]);
                            if dazi~=0
                                nn = 360/dazi + 1;
                                adp = NaN(nn,nzen+1);
                                for i=1:nn
                                    tline = fgetl(fid);
                                    linenum = linenum + 1;
                                    adp(i,:) = sscanf(tline,'%f',[1,nzen+1]);
                                end
                                atx.bds.adp{sno,fno} = adp;
                            end
                        end
                    end
                end
                
            elseif strcmp(tag,'TYPE / SERIAL NO') && strcmp(strtrim(tline(1:20)),type)
                while ~strcmp(tag,'END OF ANTENNA')
                    tline = fgetl(fid); linenum = linenum + 1;
                    tag   = strtrim(tline(61:end));
    
                    if strcmp(tag,'DAZI')
                        dazi = sscanf(tline(3:8),'%f');
                        atx.rcv.daz(1) = dazi; 
                    end
                    if strcmp(tag,'ZEN1 / ZEN2 / DZEN')
                        dzen = sscanf(tline(3:20),'%f');
                        atx.rcv.daz(2:4) = dzen';
                        nzen = (dzen(2)-dzen(1))/dzen(3) + 1;
                    end
                    if strcmp(tag,'START OF FREQUENCY') && strcmp(tline(4:6),'G01')
                        fno = 1;
                        tline = fgetl(fid);
                        tag   = strtrim(tline(61:end));
                        if strcmp(tag,'NORTH / EAST / UP')
                            atx.rcv.neu(1,:,fno) = sscanf(tline,'%f',[1,3]);
                            tline = fgetl(fid); linenum = linenum + 1;
                            atx.rcv.nad{:,fno} = sscanf(tline(9:end),'%f',[1,nzen]);
                            if dazi~=0
                                nn = 360/dazi + 1;
                                adp = NaN(nn,nzen+1);
                                for i=1:nn
                                    tline = fgetl(fid);
                                    linenum = linenum + 1;
                                    adp(i,:) = sscanf(tline,'%f',[1,nzen+1]);
                                end
                                atx.rcv.adp{:,fno} = adp;
                            end
                        end
                    elseif strcmp(tag,'START OF FREQUENCY') && strcmp(tline(4:6),'G02')
                        fno = 2;
                        tline = fgetl(fid); linenum = linenum + 1;
                        tag   = strtrim(tline(61:end));
                        if strcmp(tag,'NORTH / EAST / UP')
                            atx.rcv.neu(1,:,fno) = sscanf(tline,'%f',[1,3]);
                            tline = fgetl(fid); linenum = linenum + 1;
                            atx.rcv.nad{:,fno} = sscanf(tline(9:end),'%f',[1,nzen]);
                            if dazi~=0
                                nn = 360/dazi + 1;
                                adp = NaN(nn,nzen+1);
                                for i=1:nn
                                    tline = fgetl(fid);
                                    linenum = linenum + 1;
                                    adp(i,:) = sscanf(tline,'%f',[1,nzen+1]);
                                end
                                atx.rcv.adp{:,fno} = adp;
                            end
                        end
                    elseif strcmp(tag,'START OF FREQUENCY') && strcmp(tline(4:6),'R01')
                        fno = 3;
                        tline = fgetl(fid); linenum = linenum + 1;
                        tag   = strtrim(tline(61:end));
                        if strcmp(tag,'NORTH / EAST / UP')
                            atx.rcv.neu(1,:,fno) = sscanf(tline,'%f',[1,3]);
                            tline = fgetl(fid); linenum = linenum + 1;
                            atx.rcv.nad{:,fno} = sscanf(tline(9:end),'%f',[1,nzen]);
                            if dazi~=0
                                nn = 360/dazi + 1;
                                adp = NaN(nn,nzen+1);
                                for i=1:nn
                                    tline = fgetl(fid);
                                    linenum = linenum + 1;
                                    adp(i,:) = sscanf(tline,'%f',[1,nzen+1]);
                                end
                                atx.rcv.adp{:,fno} = adp;
                            end
                        end
                    elseif strcmp(tag,'START OF FREQUENCY') && strcmp(tline(4:6),'R02')
                        fno = 4;
                        tline = fgetl(fid); linenum = linenum + 1;
                        tag   = strtrim(tline(61:end));
                        if strcmp(tag,'NORTH / EAST / UP')
                            atx.rcv.neu(1,:,fno) = sscanf(tline,'%f',[1,3]);
                            tline = fgetl(fid); linenum = linenum + 1;
                            atx.rcv.nad{:,fno} = sscanf(tline(9:end),'%f',[1,nzen]);
                            if dazi~=0
                                nn = 360/dazi + 1;
                                adp = NaN(nn,nzen+1);
                                for i=1:nn
                                    tline = fgetl(fid);
                                    linenum = linenum + 1;
                                    adp(i,:) = sscanf(tline,'%f',[1,nzen+1]);
                                end
                                atx.rcv.adp{:,fno} = adp;
                            end
                        end
                    elseif strcmp(tag,'START OF FREQUENCY') && strcmp(tline(4:6),'E01')
                        fno = 5;
                        tline = fgetl(fid); linenum = linenum + 1;
                        tag   = strtrim(tline(61:end));
                        if strcmp(tag,'NORTH / EAST / UP')
                            atx.rcv.neu(1,:,fno) = sscanf(tline,'%f',[1,3]);
                            tline = fgetl(fid); linenum = linenum + 1;
                            atx.rcv.nad{:,fno} = sscanf(tline(9:end),'%f',[1,nzen]);
                            if dazi~=0
                                nn = 360/dazi + 1;
                                adp = NaN(nn,nzen+1);
                                for i=1:nn
                                    tline = fgetl(fid);
                                    linenum = linenum + 1;
                                    adp(i,:) = sscanf(tline,'%f',[1,nzen+1]);
                                end
                                atx.rcv.adp{:,fno} = adp;
                            end
                        end
                    elseif strcmp(tag,'START OF FREQUENCY') && strcmp(tline(4:6),'E05')
                        fno = 6;
                        tline = fgetl(fid); linenum = linenum + 1;
                        tag   = strtrim(tline(61:end));
                        if strcmp(tag,'NORTH / EAST / UP')
                            atx.rcv.neu(1,:,fno) = sscanf(tline,'%f',[1,3]);
                            tline = fgetl(fid); linenum = linenum + 1;
                            atx.rcv.nad{:,fno} = sscanf(tline(9:end),'%f',[1,nzen]);
                            if dazi~=0
                                nn = 360/dazi + 1;
                                adp = NaN(nn,nzen+1);
                                for i=1:nn
                                    tline = fgetl(fid);
                                    linenum = linenum + 1;
                                    adp(i,:) = sscanf(tline,'%f',[1,nzen+1]);
                                end
                                atx.rcv.adp{:,fno} = adp;
                            end
                        end
                    elseif strcmp(tag,'START OF FREQUENCY') && strcmp(tline(4:6),'C02')
                        fno = 7;
                        tline = fgetl(fid); linenum = linenum + 1;
                        tag   = strtrim(tline(61:end));
                        if strcmp(tag,'NORTH / EAST / UP')
                            atx.rcv.neu(1,:,fno) = sscanf(tline,'%f',[1,3]);
                            tline = fgetl(fid); linenum = linenum + 1;
                            atx.rcv.nad{:,fno} = sscanf(tline(9:end),'%f',[1,nzen]);
                            if dazi~=0
                                nn = 360/dazi + 1;
                                adp = NaN(nn,nzen+1);
                                for i=1:nn
                                    tline = fgetl(fid);
                                    linenum = linenum + 1;
                                    adp(i,:) = sscanf(tline,'%f',[1,nzen+1]);
                                end
                                atx.rcv.adp{:,fno} = adp;
                            end
                        end
                    elseif strcmp(tag,'START OF FREQUENCY') && strcmp(tline(4:6),'C06')
                        fno = 8;
                        tline = fgetl(fid); linenum = linenum + 1;
                        tag   = strtrim(tline(61:end));
                        if strcmp(tag,'NORTH / EAST / UP')
                            atx.rcv.neu(1,:,fno) = sscanf(tline,'%f',[1,3]);
                            tline = fgetl(fid); linenum = linenum + 1;
                            atx.rcv.nad{:,fno} = sscanf(tline(9:end),'%f',[1,nzen]);
                            if dazi~=0
                                nn = 360/dazi + 1;
                                adp = NaN(nn,nzen+1);
                                for i=1:nn
                                    tline = fgetl(fid);
                                    linenum = linenum + 1;
                                    adp(i,:) = sscanf(tline,'%f',[1,nzen+1]);
                                end
                                atx.rcv.adp{:,fno} = adp;
                            end
                        end
                    end
                end 
            end
        else
            continue
        end
    end

    data.atx = atx;
    fclose('all');

catch ME
    data.err.cod = 1;
    msg{1} = ['Error message:',ME.message];
    data.err.msg = msg;
    return
end
end