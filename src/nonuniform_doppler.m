function f_d = nonuniform_doppler(state, cfg)
%NONUNIFORM_DOPPLER Element-wise Doppler caused by radial/tangential motion.

channel = spherical_channel(state(1), state(2), cfg);
r = state(1);
theta = state(2);
vr = state(3);
vt = state(4);
e_r = [sin(theta); cos(theta)];
e_t = [cos(theta); -sin(theta)];
velocity = vr * e_r + vt * e_t;
delta = [channel.positions.' - channel.user_position(1); ...
    -channel.user_position(2) * ones(1, cfg.n_ant)];
unit_from_user_to_element = delta ./ channel.distances.';
range_rate = sum(unit_from_user_to_element .* (-velocity), 1).';
f_d = -(cfg.fc / cfg.c) * range_rate;
end
