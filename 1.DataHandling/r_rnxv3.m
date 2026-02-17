function [data] = r_rnxv3(fname)

data.err.cod = 0;
data.err.msg = cell(0);
head = struct;

try
    [fid,~] = fopen(fname);
    
    %%
    head.time.leap = [];

    if size(fname,2)>28
        head.time.int = str2double(fname(end-9:end-8));
    else
        head.time.int  = 1;
    end

    head.time.last = 86400;

    head.marker.name = '';
    head.rec.type = '';
    head.ant.type = '';
    head.rec.pos = NaN(1,3);
    head.time.int = NaN(1,1);

    head.seq.gps    = cell(1);
    head.seq.glo    = cell(1);
    head.seq.gal    = cell(1);
    head.seq.bds    = cell(1);
    
    head.nob.gps = 1;
    head.nob.glo = 1;
    head.nob.gal = 1;
    head.nob.bds = 1;

    while 1
        tline = fgetl(fid);
        tag  = strtrim(tline(61:end));

        switch tag
            case 'RINEX VERSION / TYPE'
                head.rinex.ver = sscanf(tline(1:20),'%f');
                head.rinex.type = sscanf(tline(21),'%c');
                if strcmp(sscanf(tline(21),'%c'),'O')
                    head.rinex.type = sscanf(tline(21),'%c');
                else
                    data.err.cod = 1;
                    data.err.msg = {'This file is not an observation file.'};
                    return
                end
                head.rinex.system = sscanf(tline(41),'%c');
    
            case 'PGM / RUN BY / DATE'
                head.agency.soft = strtrim(tline( 1:20));
                head.agency.name = strtrim(tline(21:40));
                head.agency.date = strtrim(tline(41:60));
            
            case 'MARKER NAME'
                head.marker.name = strtrim(tline(1:60));
    
            case 'OBSERVER / AGENCY'
                head.obser.name = strtrim(tline( 1:20));
                head.obser.agen = strtrim(tline(21:40));
    
            case 'REC # / TYPE / VERS'
                head.rec.number  = strtrim(tline( 1:20));
                head.rec.type = strtrim(tline(21:40));
                head.rec.version = strtrim(tline(41:60));
                
            case 'ANT # / TYPE'
                head.ant.number = strtrim(tline( 1:20));
                head.ant.type = strtrim(tline(21:40));
                
            case 'APPROX POSITION XYZ'
                head.rec.pos = sscanf(tline(1:60),'%f',[1,3]);
                
            case 'ANTENNA: DELTA H/E/N'
                head.ant.hen = sscanf(tline(1:60),'%f',[1,3]);
                
            case 'SYS / # / OBS TYPES'
                if strcmp(tline(1),'G')
                    no = sscanf(tline(5:6),'%d');
                    head.nob.gps = no;
                    if no<14
                        lst = sscanf(tline(8:60),'%s');
                    elseif no<27
                        l1 = sscanf(tline(8:60),'%s');
                        tline = fgetl(fid);
                        l2 = sscanf(tline(8:60),'%s');
                        lst = strcat(l1,l2);
                    else
                        l1 = sscanf(tline(8:60),'%s');
                        tline = fgetl(fid);
                        l2 = sscanf(tline(8:60),'%s');
                        tline = fgetl(fid);
                        l3 = sscanf(tline(8:60),'%s');
                        lst = strcat(l1,l2,l3);
                    end
                    head.seq.gps = cell(no,1);
                    for i=1:no
                        st = (i-1)*3 + 1;
                        fn = (i)*3;
                        head.seq.gps{i} = lst(st:fn);
                    end

                elseif strcmp(tline(1),'R')
                    no = sscanf(tline(5:6),'%d');
                    head.nob.glo = no;
                    if no<14
                        lst = sscanf(tline(8:60),'%s');
                    elseif no<27
                        l1 = sscanf(tline(8:60),'%s');
                        tline = fgetl(fid);
                        l2 = sscanf(tline(8:60),'%s');
                        lst = strcat(l1,l2);
                    else
                        l1 = sscanf(tline(8:60),'%s');
                        tline = fgetl(fid);
                        l2 = sscanf(tline(8:60),'%s');
                        tline = fgetl(fid);
                        l3 = sscanf(tline(8:60),'%s');
                        lst = strcat(l1,l2,l3);
                    end
                    head.seq.glo = cell(no,1);
                    for i=1:no
                        st = (i-1)*3 + 1;
                        fn = (i)*3;
                        head.seq.glo{i} = lst(st:fn);
                    end

                elseif strcmp(tline(1),'E')
                    no = sscanf(tline(5:6),'%d');
                    head.nob.gal = no;
                    if no<14
                        lst = sscanf(tline(8:60),'%s');
                    elseif no<27
                        l1 = sscanf(tline(8:60),'%s');
                        tline = fgetl(fid);
                        l2 = sscanf(tline(8:60),'%s');
                        lst = strcat(l1,l2);
                    else
                        l1 = sscanf(tline(8:60),'%s');
                        tline = fgetl(fid);
                        l2 = sscanf(tline(8:60),'%s');
                        tline = fgetl(fid);
                        l3 = sscanf(tline(8:60),'%s');
                        lst = strcat(l1,l2,l3);
                    end
                    head.seq.gal = cell(no,1);
                    for i=1:no
                        st = (i-1)*3 + 1;
                        fn = (i)*3;
                        head.seq.gal{i} = lst(st:fn);
                    end

                elseif strcmp(tline(1),'C')
                    no = sscanf(tline(5:6),'%d');
                    head.nob.bds = no;
                    if no<14
                        lst = sscanf(tline(8:60),'%s');
                    elseif no<27
                        l1 = sscanf(tline(8:60),'%s');
                        tline = fgetl(fid);
                        l2 = sscanf(tline(8:60),'%s');
                        lst = strcat(l1,l2);
                    else
                        l1 = sscanf(tline(8:60),'%s');
                        tline = fgetl(fid);
                        l2 = sscanf(tline(8:60),'%s');
                        tline = fgetl(fid);
                        l3 = sscanf(tline(8:60),'%s');
                        lst = strcat(l1,l2,l3);
                    end
                    head.seq.bds = cell(no,1);
                    for i=1:no
                        st = (i-1)*3 + 1;
                        fn = (i)*3;
                        head.seq.bds{i} = lst(st:fn);
                    end
                end
                
            case 'INTERVAL'
                head.time.int = sscanf(tline(1:10),'%f');
                
            case 'TIME OF FIRST OBS'
                head.time.first  = sscanf(tline( 1:44),'%d');
                head.time.system = sscanf(tline(45:60),'%s');
                [doy] = c_doy(head.time.first(1),head.time.first(2),head.time.first(3));
                head.time.doy = doy;
    
            case 'TIME OF LAST OBS'
                head.time.last   = sscanf(tline( 1:44),'%d');
                
            case 'LEAP SECONDS'
                head.time.leap = sscanf(tline(1:6),'%d');
                
            case 'END OF HEADER'
                break
        end
    end

    % if isempty(head.ant.type)
    %     data.err.cod = 1;
    %     data.err.msg = {'The antenna name does not exist.'};
    %     return
    % end
    % if isnan(head.rec.pos)
    %     data.err.cod = 1;
    %     data.err.msg = {'The approximate position does not exist.'};
    %     return
    % end
    % if isnan(head.time.int)
    %     data.err.cod = 1;
    %     data.err.msg = {'The observation interval does not exist.'};
    %     return
    % end

    if isempty(head.time.leap)
        [~,mjd] = cal2jul(head.time.first(1),head.time.first(2),head.time.first(3),...
            (head.time.first(4)*3600 + head.time.first(5)*60 + head.time.first(6)));
        
        head.time.leap = d_leap(mjd);
        if strcmp(head.time.system,'GPS')
            head.time.leap = head.time.leap - 19;
        end
    end

    data.rnx.head = head;
    
    %% 
    if head.time.last == 86400
        fi = head.time.first(4,1)*3600 + head.time.first(5,1)*60 + head.time.first(6,1);
        la = head.time.last;
        max = round((la-fi)/head.time.int + 1);
    else
        fi = head.time.first(4,1)*3600 + head.time.first(5,1)*60 + head.time.first(6,1);
        la = head.time.last (4,1)*3600 + head.time.last (5,1)*60 + head.time.last (6,1);
        max = round((la-fi)/head.time.int + 1);
    end

    gobs = NaN(max,32,head.nob.gps); 
    robs = NaN(max,26,head.nob.glo);
    eobs = NaN(max,36,head.nob.gal);
    cobs = NaN(max,62,head.nob.bds);

    eps  = NaN(max, 1);

    epno = 0;
    lno  = 0;
    while ~feof(fid)
        tline = fgetl(fid);
        lno   = lno + 1;

        if ~isempty(tline) && strcmp(tline(1),'>')
            nep  = sscanf(tline(3:end),'%f');
            epno = epno + 1;
            eps(epno,1) = nep(4)*3600 + nep(5)*60 + nep(6); 
            for i=1:nep(8)
                tline = fgetl(fid);
                lno   = lno + 1;

                if strcmp(tline(1),'G')
                    k   = sscanf(tline(2:3),'%d'); 
                    if k>size(gobs,2)
                        continue
                    else
                        nu  = 0;
                        for u=4:16:size(tline,2)
                            nu = nu + 1;
                            tls = sscanf(tline(u:u+13),'%f');
                            if numel(tls)>1
                                break
                            elseif ~isempty(tls)
                                if strcmp(head.seq.gps{nu}(1),'L')
                                    ff = str2double(head.seq.gps{nu}(2));
                                    [~,wavl] = gnss_freq('G',k,ff);
                                    gobs(epno,k,nu) = tls*wavl;
                                else
                                    gobs(epno,k,nu) = tls;
                                end
                            end
                        end
                    end

                elseif strcmp(tline(1),'R')
                    k   = sscanf(tline(2:3),'%d');
                    if k>24
                        continue
                    else
                        nu  = 0;
                        for u=4:16:size(tline,2)
                            nu = nu + 1;
                            tls = sscanf(tline(u:u+13),'%f');
                            if numel(tls)>1
                                break
                            elseif ~isempty(tls)
                                if strcmp(head.seq.glo{nu}(1),'L')
                                    ff = str2double(head.seq.glo{nu}(2));
                                    [~,wavl] = gnss_freq('R',k,ff);
                                    robs(epno,k,nu) = tls*wavl;
                                else
                                    robs(epno,k,nu) = tls;
                                end
                            end
                        end
                    end

                elseif strcmp(tline(1),'E')
                    k   = sscanf(tline(2:3),'%d');
                    if k>size(eobs,2)
                        continue
                    else
                        nu  = 0;
                        for u=4:16:size(tline,2)
                            nu = nu + 1;
                            tls = sscanf(tline(u:u+13),'%f');
                            if numel(tls)>1
                                break
                            elseif ~isempty(tls)
                                if strcmp(head.seq.gal{nu}(1),'L')
                                    ff = str2double(head.seq.gal{nu}(2));
                                    [~,wavl] = gnss_freq('E',k,ff);
                                    eobs(epno,k,nu) = tls*wavl;
                                else
                                    eobs(epno,k,nu) = tls;
                                end
                            end
                        end
                    end

                elseif strcmp(tline(1),'C')
                    k   = sscanf(tline(2:3),'%d');
                    if k>size(cobs,2)
                        continue
                    else
                        nu  = 0;
                        for u=4:16:size(tline,2)
                            nu = nu + 1;
                            tls = sscanf(tline(u:u+13),'%f');
                            if numel(tls)>1
                                break
                            elseif ~isempty(tls)
                                if strcmp(head.seq.bds{nu}(1),'L')
                                    ff = str2double(head.seq.bds{nu}(2));
                                    [~,wavl] = gnss_freq('C',k,ff);
                                    cobs(epno,k,nu) = tls*wavl;
                                else
                                    cobs(epno,k,nu) = tls;
                                end
                            end
                        end
                    end

                end
            end
        end
    end

    if max>epno
        gobs(epno+1:max,:,:) = [];
        robs(epno+1:max,:,:) = [];
        eobs(epno+1:max,:,:) = [];
        cobs(epno+1:max,:,:) = [];
        eps(epno+1:max,:) = [];
    end

    data.rnx.obs.gps = gobs;
    data.rnx.obs.glo = robs;
    data.rnx.obs.gal = eobs;
    data.rnx.obs.bds = cobs;
    data.rnx.obs.ep  = eps;

    fclose('all');
catch ME
    data.err.cod = 1;
    msg{1} = ['Error message:',ME.message];
    data.err.msg = msg;
    return
end

end