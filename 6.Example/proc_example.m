files.rnx = 'FLY_0098.obs';
files.sp3 = 'COD0MGXFIN_20223090000_01D_05M_ORB.SP3';
files.clk = 'COD0MGXFIN_20223090000_01D_30S_CLK.CLK';
files.atx = 'igs20_2233.atx';
files.osb = 'COD0MGXFIN_20223090000_01D_01D_OSB.BIA';

options.clkint = 30;
options.elvang = 10;

options.SatClk = 1;
options.SatAPC = 1;
options.RecAPC = 1;
options.RecARP = 1;
options.RelClk = 1;
options.SatWind = 1;
options.AtmTrop = 1;
options.RelPath = 1;
options.Solid = 1;

options.vrt.gps.p = 1; options.vrt.gps.l = 1;
options.vrt.glo.p = 2; options.vrt.glo.l = 1;
options.vrt.gal.p = 1; options.vrt.gal.l = 1;
options.vrt.bds.p = 1; options.vrt.bds.l = 1;

options.ProMod = 1;

options.IntPos   =  0;
options.IntPos2  =  0;
options.IntClk   =  1;
options.IntClk2  =  5;
options.IntTrop  =  0.5;
options.IntTrop2 =  0;
options.IntAmb   =  2;
options.IntAmb2  =  1;

options.NosPos   =  1;
options.NosPos2  = -1;
options.NosClk   =  1;
options.NosClk2  =  5;
options.NosTrop  =  1;
options.NosTrop2 = -7;

options.CodeStd  = 3;
options.PhaseStd = 0.003;

[data,options] = datahandling(files,options);

[data] = preprocessing(data,options);

[xs,pks,kofs,ress,satno,fxpos,tfxl] = modelandfiltering(data,options);