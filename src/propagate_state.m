function next = propagate_state(state, cfg, add_noise, include_hardware)
%PROPAGATE_STATE Cartesian constant-velocity motion with polar state output.
% The tracker state is [r, theta, v_r, v_t, CFO, CPE], but the physical
% motion is advanced in Cartesian coordinates.  This avoids singular-looking
% polar dynamics and unrealistic centripetal terms when the user is close to
% the array in a strong near-field setting.

if nargin < 3
    add_noise = false;
end
if nargin < 4
    include_hardware = true;
end

dt = cfg.track_period;
r = max(state(1), cfg.channel.min_range);
theta = state(2);
vr = state(3);
vt = state(4);

next = state;
e_r = [sin(theta); cos(theta)];
e_t = [cos(theta); -sin(theta)];
position = r * e_r;
velocity = vr * e_r + vt * e_t;

if add_noise
    local_acceleration = cfg.motion.accel_std * randn(2, 1);
    velocity = velocity + (local_acceleration(1) * e_r + local_acceleration(2) * e_t) * dt;
end

position = position + velocity * dt;
new_r = max(norm(position), cfg.channel.min_range);
if norm(position) < cfg.channel.min_range
    position = cfg.channel.min_range * e_r;
end
new_theta = atan2(position(1), position(2));
new_e_r = [sin(new_theta); cos(new_theta)];
new_e_t = [cos(new_theta); -sin(new_theta)];

next(1) = new_r;
next(2) = wrap_angle(new_theta);
next(3) = dot(velocity, new_e_r);
next(4) = dot(velocity, new_e_t);

if include_hardware
    next(6) = wrap_angle(state(6) + 2 * pi * state(5) * dt);
else
    next(5:6) = 0;
end

if add_noise
    if include_hardware
        next(5) = next(5) + cfg.motion.cfo_rw_std * randn;
        phase_std = sqrt(2 * pi * cfg.hardware.phase_noise_linewidth_hz * dt);
        next(6) = wrap_angle(next(6) + phase_std * randn);
    end
end
end
