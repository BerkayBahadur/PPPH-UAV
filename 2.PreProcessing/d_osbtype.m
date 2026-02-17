function [obb] = d_osbtype(data,sys,soss)

obb = struct;

if strcmp(sys,'G')
    % P1
    switch soss{1}
        case 'C1C'
            obb.ob1 = data.osb.gps.c1c;
        case 'C1L'
            obb.ob1 = data.osb.gps.c1l;
        case 'C1X'
            obb.ob1 = data.osb.gps.c1x;
        case 'C1W'
            obb.ob1 = data.osb.gps.c1w;
    end
    % L1
    switch soss{2}
        case 'L1C'
            obb.ob2 = data.osb.gps.l1c;
        case 'L1L'
            obb.ob2 = data.osb.gps.l1l;
        case 'L1X'
            obb.ob2 = data.osb.gps.l1x;
        case 'L1W'
            obb.ob2 = data.osb.gps.l1w;
    end
    % P2
    switch soss{3}
        case 'C2L'
            obb.ob3 = data.osb.gps.c2l;
        case 'C2S'
            obb.ob3 = data.osb.gps.c2s;
        case 'C2X'
            obb.ob3 = data.osb.gps.c2x;
        case 'C2W'
            obb.ob3 = data.osb.gps.c2w;
    end
    % L2
    switch soss{4}
        case 'L2L'
            obb.ob4 = data.osb.gps.l2l;
        case 'L2S'
            obb.ob4 = data.osb.gps.l2s;
        case 'L2X'
            obb.ob4 = data.osb.gps.l2x;
        case 'L2W'
            obb.ob4 = data.osb.gps.l2w;
    end

elseif strcmp(sys,'R')
    % P1
    switch soss{1}
        case 'C1P'
            obb.ob1 = data.osb.glo.c1p;
        case 'C1C'
            obb.ob1 = data.osb.glo.c1p;
    end
    % L1
    switch soss{2}
        case 'L1P'
            obb.ob2 = data.osb.glo.l1p;
        case 'L1C'
            obb.ob2 = data.osb.glo.l1p;
    end
    % P2
    switch soss{3}
        case 'C2P'
            obb.ob3 = data.osb.glo.c2p;
        case 'C2C'
            obb.ob3 = data.osb.glo.c2p;
    end
    % L2
    switch soss{4}
        case 'L2P'
            obb.ob4 = data.osb.glo.l2p;
        case 'L2C'
            obb.ob4 = data.osb.glo.l2p;
    end

elseif strcmp(sys,'E')
    % P1
    switch soss{1}
        case 'C1C'
            obb.ob1 = data.osb.gal.c1c;
        case 'C1X'
            obb.ob1 = data.osb.gal.c1x;
    end
    % L1
    switch soss{2}
        case 'L1C'
            obb.ob2 = data.osb.gal.l1c;
        case 'L1X'
            obb.ob2 = data.osb.gal.l1x;
    end
    % P2
    switch soss{3}
        case 'C5I'
            obb.ob3 = data.osb.gal.c5i;
        case 'C5Q'
            obb.ob3 = data.osb.gal.c5q;
        case 'C5X'
            obb.ob3 = data.osb.gal.c5x;
    end
    % L2
    switch soss{4}
        case 'L5I'
            obb.ob4 = data.osb.gal.l5i;
        case 'L5Q'
            obb.ob4 = data.osb.gal.l5q;
        case 'L5X'
            obb.ob4 = data.osb.gal.l5x;
    end

elseif strcmp(sys,'C')
    % P1
    switch soss{1}
        case 'C2I'
            obb.ob1 = data.osb.bds.c2i;
        case 'C2Q'
            obb.ob1 = data.osb.bds.c2i;
        case 'C2X'
            obb.ob1 = data.osb.bds.c2i;
    end
    % L1
    switch soss{2}
        case 'L2I'
            obb.ob2 = data.osb.bds.l2i;
        case 'L2Q'
            obb.ob2 = data.osb.bds.l2i;
        case 'L2X'
            obb.ob2 = data.osb.bds.l2i;
    end
    % P2
    switch soss{3}
        case 'C6I'
            obb.ob3 = data.osb.bds.c6i;
        case 'C6Q'
            obb.ob3 = data.osb.bds.c6i;
        case 'C6X'
            obb.ob3 = data.osb.bds.c6i;
    end
    % L2
    switch soss{4}
        case 'L6I'
            obb.ob4 = data.osb.bds.l6i;
        case 'L6Q'
            obb.ob4 = data.osb.bds.l6i;
        case 'L6X'
            obb.ob4 = data.osb.bds.l6i;
    end
end


end