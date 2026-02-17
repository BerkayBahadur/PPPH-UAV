function [arc] = d_arcs(data)

sn = size(data.rnx.pob,2);
arc = cell(1,sn);

cc = 30;

for k=1:sn
    ls = data.rnx.pob(:,k);
    es = data.rnx.obs.ep;
    % rounding depending on the interval
    N = num2str(data.rnx.head.time.int);
    rr = numel(N) - find(N == '.');
    if isempty(rr)
        ns = find(round(diff(es,1))>300); % five minutes
    else
        ns = find(round(diff(es,1),rr)>300);
    end
    if ~isempty(ns)
        ls(1+ns,1) = 0;
    end

    row = find(ls==1);
    brk = find(diff(row)~= 1);

    frs = [1;(brk)+1]; lst = [brk;size(row,1)];
    als = [frs lst];

    if size((als),1)>1
        for t=size((als),1):-1:1
            if (als(t,2)-als(t,1))<cc
                als(t,:) = [];
            end
        end
    else
        if (als(1,2)-als(1,1))<cc
            als(1,:) = [];
        end
    end

    arn = NaN(size(als,1),2);
    for i=1:size(als,1)
        arn(i,1) = row(als(i,1)); arn(i,2) = row(als(i,2));
    end

    arc{k} = arn;
end

end