function [data] = d_cycslp(data)

[arc] = d_arcs(data);
sn = size(data.rnx.pob,2);

c   = 299792458;% m/s
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

if int<=1
    bmw = 2.5;
elseif int<=20
    bmw  = 2.5/20*int + 2.5;
elseif int<=60
    bmw = 5.0;
else
    bmw  = 7.5;
end

[frqs,~] = frequencies;

for k=1:sn
    karc = arc{k};
    
    gfc = data.rnx.sob(:,k,2) - data.rnx.sob(:,k,4);% in meters
    
    f1 = frqs(k,1); f2 = frqs(k,2);
    lwl = (data.rnx.sob(:,k,2).*f1 - data.rnx.sob(:,k,4).*f2)./(f1-f2);
    pnl = (data.rnx.sob(:,k,1).*f1 + data.rnx.sob(:,k,3).*f2)./(f1+f2);
    lamwl = c/(f1-f2);
    mwc = (lwl - pnl)./(lamwl);% in cycles

    for s=1:size(karc,1)
        st = karc(s,1); fn = karc(s,2);

        for i=st+1:fn
            elk = data.sat.lcl(i,k,2);
            cgf = gfc(i,1) - gfc(i-1,1);

            % threshold for GF combination
            if elk>15
                tgf = bgf;
            else
                tgf = (-1/15*elk + 2)*bgf;
            end
            % check for GFC
            if abs(cgf)>tgf
                sgf = 1;
            else
                sgf = 0;
            end

            %
            cmw = mwc(i,1) - mwc(i-1,1); % in cycle
            
            % threshold for MW combination
            if elk>20
                tmw = bmw;
            else
                tmw= (-0.1*elk + 3)*bmw;
            end
            % check for MWC
            if abs(cmw)>tmw
                smw = 1;
            else
                smw = 0;
            end

            %
            if sgf==1 || smw==1
                data.rnx.pob(i,k) = 0;
            end
        end
    end
end

end