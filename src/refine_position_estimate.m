function [position, metric] = refine_position_estimate(y, reference_state, cfg)
%REFINE_POSITION_ESTIMATE Narrow local range-angle beam scan for tracking.
% This is used after the EKF correction in strong near-field/narrow-beam
% settings.  It models a small local refinement around the predicted focus.

snapshot = y(:, 1);
snapshot_energy = max(real(snapshot' * snapshot), eps);
ranges = max(cfg.channel.min_range, reference_state(1) + ...
    (-cfg.tracker.refine_range_half_width_m:cfg.tracker.refine_range_step_m:cfg.tracker.refine_range_half_width_m));
angles = wrap_angle(reference_state(2) + ...
    (-cfg.tracker.refine_angle_half_width_rad:cfg.tracker.refine_angle_step_rad:cfg.tracker.refine_angle_half_width_rad));

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
