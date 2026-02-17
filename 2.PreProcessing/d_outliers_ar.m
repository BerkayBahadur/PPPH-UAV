function [data] = d_outliers_ar(data)

[arc] = d_arcs(data);
sn = size(data.rnx.pob,2);

c  = 299792458; % m/s
[frqs,~] = frequencies;

for k=1:sn
    % observations
    f1 = frqs(k,1); f2 = frqs(k,2);
    pnl = (data.rnx.sob(:,k,1).*f1 + data.rnx.sob(:,k,3).*f2)./(f1+f2);
    lwl = (data.rnx.sob(:,k,2).*f1 - data.rnx.sob(:,k,4).*f2)./(f1-f2);
    %mw = (lwl - pnl);
    lamwl = c/(f1-f2);
    mw = (lwl - pnl)./(lamwl);
    % arcs
    ark = arc{k};
    for n=1:size(ark,1)
        st = ark(n,1); fn = ark(n,2);

        % 
        al = mw(st:fn);
        mm = mean(al,'omitmissing');
        %rm = round(mm);
        dm = abs(al - mm);
        rt = sum(dm>0.33)/length(dm);
        if rt>0.3
            data.rnx.pob(st:fn,k) = 0;
        end
    end
end
end