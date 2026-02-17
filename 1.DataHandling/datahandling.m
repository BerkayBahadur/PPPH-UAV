function [data,options] = datahandling(files,options)

[data] = r_rnx(files.rnx);

[data] = r_sp3(files.sp3,data);

[data] = r_clk(files.clk,data,options.clkint);

[data] = r_atx(files.atx,data);

[data] = r_osb(files.osb,data);

end