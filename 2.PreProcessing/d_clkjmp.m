function [data] = d_clkjmp(data)

[arc] = d_arcs(data);
en = size(data.rnx.pob,1);
sn = size(data.rnx.pob,2);

[frqs,~] = frequencies;

c   = 299792458;% m/s
pot  = NaN(en,sn);
val  = NaN(en,sn);
jump = zeros(en,1);

int = data.rnx.head.time.int;
if int<=1
    bgf = 0.05;
elseif int<=20
    bgf = 0.1/20*int + 0.05;
elseif int<=60
    bgf = 0.15;
elseif int<=100
    bgf = 0.25;
else
    bgf = 0.35;
end

for k=1:sn
    % frequencies and coefficients
    f1 = frqs(k,1); f2 = frqs(k,2);
    c1 = (f1^2)/(f1^2-f2^2); c2 = (f2^2)/(f1^2-f2^2);
    % observations
    ifp = c1.*data.rnx.sob(:,k,1) - c2*data.rnx.sob(:,k,3);
    ifl = c1.*data.rnx.sob(:,k,2) - c2*data.rnx.sob(:,k,4);
    dfi = ifp - ifl;
    gfl = data.rnx.sob(:,k,2) - data.rnx.sob(:,k,4);
    % arcs
    ark = arc{k};
    for t=1:size(ark,1)
        st = ark(t,1); 
        fn = ark(t,2);
        
        % epochs
        for ii=st+1:fn
            % elevation angle for the related satellite
            elk = data.sat.lcl(ii,k,2);
            % threshold for GF combination
            if elk>15
                tgf = bgf;
            else
                tgf = (-1/15*elk + 2)*bgf;
            end
            dgf = gfl(ii) - gfl(ii-1);
            % check for GFC
            if abs(dgf)>tgf
                sgf = 1;
            else
                sgf = 0;
            end
            % if no cycle slip exist
            if sgf==0
                % standard deviation for the related satellite
                thr = (10^-3*c) - (3*5);
                ddf = dfi(ii) - dfi(ii-1);
                val(ii,k) = ddf;
                if abs(ddf)>=thr
                    pot(ii,k) = 1;
                else
                    pot(ii,k) = 0;
                end
            end
        end
    end   
end

for i=1:en
    if any(~isnan(pot(i,:)))
        na = sum(~isnan(pot(i,:)));
        nn = sum(pot(i,:)==1);
    
        if na==nn
            M = mean(val(i,pot(i,:)==1),'omitmissing')*(10^3/c);

            if abs(M - round(M)) <= (10^-5) %microsecond
                jump(i) = round(M);
            end
        end
    end
end


if any(jump~=0)
    kern2 = find(jump~=0);
    for k=1:sn
        ark = arc{k};
        for t=1:size(ark,1)
            st = ark(t,1); 
            fn = ark(t,2);
            for ke = kern2'
                if (ke>st) && (ke<=fn)
                    data.rnx.sob(ke:fn,k,2) = data.rnx.sob(ke:fn,k,2) + jump(ke,1)*(c*10^-3);
                    data.rnx.sob(ke:fn,k,4) = data.rnx.sob(ke:fn,k,4) + jump(ke,1)*(c*10^-3);
                end
            end
        end
    end
end

end