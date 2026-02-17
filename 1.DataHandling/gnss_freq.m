function [freq,wavl] = gnss_freq(sys,prn,fno)

c = 299792458;

glok = [1 -4 5 6 1 -4 5 6 -2 -7 0 -1 -2 -7 0 -1 4 -3 3 2 4 -3 3 2 0 0];

switch sys
    case 'G'
        if fno==1
            freq = 10.23*10^6*154;
            wavl = c/freq;
        elseif fno==2
            freq = 10.23*10^6*120;
            wavl = c/freq;
        elseif fno==5
            freq = 10.23*10^6*115;
            wavl = c/freq;
        end
    case 'R'
        if fno==1
            k = 9/16;
            freq = (1602 + k*glok(prn))*10^6;
            wavl = c/freq;
        elseif fno==2
            k = 7/16;
            freq = (1246 + k*glok(prn))*10^6;
            wavl = c/freq;
        elseif fno==3
            freq = 10.23*10^6*117.5;
            wavl = c/freq;
        elseif fno==4
            freq = 10.23*10^6*156.5;
            wavl = c/freq;
        elseif fno==6
            freq = 10.23*10^6*122;
            wavl = c/freq;
        end
    case 'E'
        if fno==1
            freq = 10.23*10^6*154;
            wavl = c/freq;
        elseif fno==5
            freq = 10.23*10^6*115;
            wavl = c/freq;
        elseif fno==6
            freq = 10.23*10^6*125;
            wavl = c/freq;
        elseif fno==7
            freq = 10.23*10^6*118;
            wavl = c/freq;
        elseif fno==8
            freq = 10.23*10^6*116.5;
            wavl = c/freq;
        end
    case 'C'
        if fno==1
            freq = 10.23*10^6*154;
            wavl = c/freq;
        elseif fno==2
            freq = 10.23*10^6*152.6;
            wavl = c/freq;
        elseif fno==5
            freq = 10.23*10^6*115;
            wavl = c/freq;
        elseif fno==6
            freq = 10.23*10^6*124;
            wavl = c/freq;
        elseif fno==7
            freq = 10.23*10^6*118;
            wavl = c/freq;
        elseif fno==8
            freq = 10.23*10^6*116.5;
            wavl = c/freq;
        end
    case 'S'
        if fno==1
            freq = 10.23*10^6*154;
            wavl = c/freq;
        elseif fno==5
            freq = 10.23*10^6*115;
            wavl = c/freq;
        end
end