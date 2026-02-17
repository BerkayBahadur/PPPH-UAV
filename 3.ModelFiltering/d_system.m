function [bp,ap,sit] = d_system(data)

all = [1 0 0 0;%G
       0 1 0 0;%R
       0 0 1 0;%E
       0 0 0 1;%C
       1 1 0 0;%GR
       1 0 1 0;%GE
       1 0 0 1;%GC
       0 1 1 0;%RE
       0 1 0 1;%RC
       0 0 1 1;%EC
       1 1 1 0;%GRE
       1 1 0 1;%GRC
       1 0 1 1;%GEC
       0 1 1 1;%REC
       1 1 1 1;];%GREC

sns = [32 26 36 62];

pss = [ 1  1  1  2];

cmp = [0 0 0 0];
cmp(1) = any(any(data.rnx.pob(:,1:32)==1));
cmp(2) = any(any(data.rnx.pob(:,33:58)==1));
cmp(3) = any(any(data.rnx.pob(:,59:94)==1));
cmp(4) = any(any(data.rnx.pob(:,95:156)==1));

for i=1:size(all,1)
    if isequal(cmp,all(i,:))
        sit = i;
        ap = cmp*sns';
        bp = 4 + cmp*pss';
        break
    end
end

end