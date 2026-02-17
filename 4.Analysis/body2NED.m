function offsetNED = body2NED(yaw, pitch, roll, offsetBody)

yaw = deg2rad(yaw);
pitch = deg2rad(pitch);
roll = deg2rad(roll);

Rz_yaw = [
    cos(yaw), -sin(yaw), 0;
    sin(yaw),  cos(yaw), 0;
    0,         0,        1
    ];

Ry_pitch = [
    cos(pitch),  0, sin(pitch);
    0,           1, 0;
    -sin(pitch), 0, cos(pitch)
    ];

Rx_roll = [
    1, 0,         0;
    0, cos(roll), -sin(roll);
    0, sin(roll),  cos(roll)
    ];


Cn_b = Rz_yaw * Ry_pitch * Rx_roll;


offsetNED = (Cn_b * offsetBody')';

end