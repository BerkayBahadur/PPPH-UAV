function [data] = preprocessing(data,options)

[data] = s_obs(data);

[data] = c_satpos(data);

[data] = c_local(data,options.elvang);

[data] = d_clkjmp(data);

[data] = d_outliers(data);

[data] = d_cycslp(data);

%[data] = c_osb(data);

[data] = c_arcobs(data);

[data] = d_outliers_ar(data);

end