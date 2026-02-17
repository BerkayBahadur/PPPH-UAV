function [freq,wavl] = frequencies

sno = 156;
freq = zeros(sno,2);
wavl = zeros(sno,2);


for k=1:sno
    if k<33
        prn = k;
        [f1,w1] = gnss_freq('G',prn,1);
        [f2,w2] = gnss_freq('G',prn,2);
    elseif k<59
        prn = k - 32;
        [f1,w1] = gnss_freq('R',prn,1);
        [f2,w2] = gnss_freq('R',prn,2);
    elseif k<95
        prn = k - 58;
        [f1,w1] = gnss_freq('E',prn,1);
        [f2,w2] = gnss_freq('E',prn,5);
    else
        prn = k - 94;
        [f1,w1] = gnss_freq('C',prn,2);
        [f2,w2] = gnss_freq('C',prn,6);
    end
    freq(k,:) = [f1 f2];
    wavl(k,:) = [w1 w2];
end

end