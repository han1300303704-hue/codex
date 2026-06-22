function [position, metric] = coarse_position_estimate(y, reference_state, cfg)
%COARSE_POSITION_ESTIMATE Pilot-aided local range-angle acquisition.
% A common CFO/CPE rotates all array elements equally, so a coherent
% beam-scan magnitude on the first pilot is suitable before EKF start-up.

snapshot = y(:, 1);
snapshot_energy = max(real(snapshot' * snapshot), eps);
ranges = max(cfg.channel.min_range, reference_state(1) + ...
    (-cfg.initializer.range_half_width_m:cfg.initializer.range_step_m:cfg.initializer.range_half_width_m));
angles = wrap_angle(reference_state(2) + ...
    (-cfg.initializer.angle_half_width_rad:cfg.initializer.angle_step_rad:cfg.initializer.angle_half_width_rad));

metric = -inf;
position = [reference_state(1); reference_state(2)];
for r = ranges
    for theta = angles
        candidate = focus_beam(r, theta, inf, cfg);
        score = abs(candidate.weights' * snapshot)^2 / snapshot_energy;
        if score > metric
            metric = score;
            position = [r; theta];
        end
    end
end
end
