function [data] = r_rnx(fname)

[ver] = get_ver(fname);

if floor(ver)==3
    [data] = r_rnxv3(fname);
else
    data.err.cod = 1;
    data.err.msg = {'RINEX version should be 3.xx.'};
    return
end

end

function [ver] = get_ver(fname)

[fid,errmsg] = fopen(fname);

if ~any(errmsg)
    while 1
        tline = fgetl(fid);
        tag  = strtrim(tline(61:end));
        if strcmp(tag,'RINEX VERSION / TYPE')
            ver = sscanf(tline(1:20),'%f');
            break
        end
    end
end

fclose('all');
end