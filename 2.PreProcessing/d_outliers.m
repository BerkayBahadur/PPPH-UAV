function [data] = d_outliers(data)

[arc] = d_arcs(data);
sn = size(data.rnx.pob,2);

%c  = 299792458; % m/s
[frqs,~] = frequencies;

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
    % observations
    f1 = frqs(k,1); f2 = frqs(k,2);
    pnl = (data.rnx.sob(:,k,1).*f1 + data.rnx.sob(:,k,3).*f2)./(f1+f2);
    lwl = (data.rnx.sob(:,k,2).*f1 - data.rnx.sob(:,k,4).*f2)./(f1-f2);
    mw = (lwl - pnl);
    %lamwl = c/(f1-f2);
    % nwl = (lwl - pnl)./(lamwl);
    gfl = data.rnx.sob(:,k,2) - data.rnx.sob(:,k,4);
    % arcs
    ark = arc{k};
    for n=1:size(ark,1)
        st = ark(n,1); fn = ark(n,2);
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
                %smw = std(mw(st:ii),'omitmissing');
                % measure
                dmw = mw(ii) - mw(ii-1);
                % chechk
                if abs(dmw)>=2%(5*smw)
                    data.rnx.pob(ii,k) = 0;
                    %mw(data.rnx.pob(:,k)==0) = NaN;
                end
            end
        end
    end
end
end