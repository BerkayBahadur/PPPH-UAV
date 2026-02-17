function pppar_sol = process_PPPAR_solutions(input_filename, data, xs, offsetBody, output_filename)

data_table = readtable(input_filename);

ypr_data.image_number = data_table.image_name;
ypr_data.yaw = data_table.yaw;
ypr_data.pitch = data_table.pitch;
ypr_data.roll = data_table.roll;


pppar_pos.gps_timestamp = [];
pppar_pos.gps_week = [];
[c_lat,c_lon,pppar_pos.h] = xyz2ell(xs(:,1),xs(:,2),xs(:,3));
pppar_pos.lat = rad2deg(c_lat);
pppar_pos.lon = rad2deg(c_lon);


fractional_hour = data.rnx.obs.ep/3600;
day_of_year = data.rnx.head.time.doy*ones(size(data.rnx.obs.ep));
year = data.rnx.head.time.first(1)*ones(size(data.rnx.obs.ep));


for i = 1:length(year)
    gps_time = datetime(year(i), 1, 1) + day_of_year(i) - 1 + fractional_hour(i) / 24;
    gps_epoch = datetime(1980, 1, 6);
    seconds_since_epoch = (gps_time - gps_epoch) * 24 * 3600;

    gps_week = floor(seconds_since_epoch / (7 * 24 * 3600));
    gps_seconds_of_week = mod(seconds_since_epoch, 7 * 24 * 3600);

    pppar_pos.gps_timestamp = [pppar_pos.gps_timestamp; gps_seconds_of_week];
    pppar_pos.gps_week = [pppar_pos.gps_week; gps_week];
end


pppar_sol.gps_timestamp = data_table.gps_timestamp;
pppar_sol.gps_week = data_table.gps_week;
pppar_sol.lat_uncorrected = interp1(pppar_pos.gps_timestamp, pppar_pos.lat, data_table.gps_timestamp, 'spline', 'extrap');
pppar_sol.lon_uncorrected = interp1(pppar_pos.gps_timestamp, pppar_pos.lon, data_table.gps_timestamp, 'spline', 'extrap');
pppar_sol.h_uncorrected = interp1(pppar_pos.gps_timestamp, pppar_pos.h, data_table.gps_timestamp, 'spline', 'extrap');


offsetNED = NaN(size(ypr_data.yaw,1),3);
for i = 1:size(offsetNED,1)
    offsetNED(i,:) = body2NED(ypr_data.yaw(i), ypr_data.pitch(i), ypr_data.roll(i), offsetBody);
end
offsetNEU = [offsetNED(:,1:2) -offsetNED(:,3)];


offsetXYZ = NaN(size(offsetNEU,1),3);
for j = 1:size(offsetXYZ,1)
    offsetXYZ(j,:) = NEU2xyz(offsetNEU(j,:),deg2rad(pppar_sol.lat_uncorrected(j)),deg2rad(pppar_sol.lon_uncorrected(j)));
end


[pppar_sol.x_uncorrected, pppar_sol.y_uncorrected, pppar_sol.z_uncorrected] = plh2xyz(pppar_sol.lat_uncorrected, pppar_sol.lon_uncorrected, pppar_sol.h_uncorrected);


pppar_sol.x_corrected = pppar_sol.x_uncorrected + offsetXYZ(:,1);
pppar_sol.y_corrected = pppar_sol.y_uncorrected + offsetXYZ(:,2);
pppar_sol.z_corrected = pppar_sol.z_uncorrected + offsetXYZ(:,3);

[lat_c,lon_c,pppar_sol.h_corrected] = xyz2ell(pppar_sol.x_corrected,pppar_sol.y_corrected,pppar_sol.z_corrected);
pppar_sol.lat_corrected = rad2deg(lat_c);
pppar_sol.lon_corrected = rad2deg(lon_c);

if nargin == 5 && ~isempty(output_filename)
    fileID = fopen(output_filename, 'w');
    
    for i = 1:length(ypr_data.image_number)
        fprintf(fileID, '%s,%.8f,%.8f,%.3f,%.1f,%.1f,%.1f\n', ...
            ypr_data.image_number{i}, ...
            pppar_sol.lon_corrected(i), ...
            pppar_sol.lat_corrected(i), ...
            pppar_sol.h_corrected(i), ...
            0, 0, 0);
    end

    fclose(fileID);
end
end