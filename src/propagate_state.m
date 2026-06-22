function next = propagate_state(state, cfg, add_noise, include_hardware)
%PROPAGATE_STATE Polar free-motion model with CFO and phase evolution.

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
next(1) = max(cfg.channel.min_range, r + vr * dt);
next(2) = wrap_angle(theta + (vt / r) * dt);
next(3) = vr + (vt^2 / r) * dt;
next(4) = vt - (vr * vt / r) * dt;

if include_hardware
    next(6) = wrap_angle(state(6) + 2 * pi * state(5) * dt);
else
    next(5:6) = 0;
end

if add_noise
    acceleration = cfg.motion.accel_std * randn(2, 1);
    next(3:4) = next(3:4) + acceleration * dt;
    if include_hardware
        next(5) = next(5) + cfg.motion.cfo_rw_std * randn;
        phase_std = sqrt(2 * pi * cfg.hardware.phase_noise_linewidth_hz * dt);
        next(6) = wrap_angle(next(6) + phase_std * randn);
    end
end
end
