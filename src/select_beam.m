function [beam, selection] = select_beam(x_hat, P, y, offsets, cfg, use_quantization, robust_selection, force_reacquire)
%SELECT_BEAM Choose a continuous/quantized focusing beam for the next slot.
% Robust local selection maximizes uncertainty-averaged gain and penalizes
% guard-ring response, which suppresses quantization-induced pseudo-focuses.

if nargin < 8
    force_reacquire = false;
end
if use_quantization
    bits = cfg.hardware.phase_shifter_bits;
else
    bits = inf;
end

if force_reacquire
    [beam, selection] = reacquire_beam(y, cfg, bits);
    return;
end

if ~robust_selection
    beam = focus_beam(x_hat(1), x_hat(2), bits, cfg);
    selection = struct('mode', 'predicted_center', 'score', NaN, ...
        'candidate', beam.focus, 'num_candidates', 1);
    return;
end

sigma_r = max([sqrt(P(1, 1)), cfg.beam.robust_min_range_std_m, 0.02]);
sigma_theta = max([sqrt(P(2, 2)), cfg.beam.robust_min_angle_std_rad, deg2rad(0.02)]);
offsets_local = cfg.beam.local_offsets * cfg.beam.local_sigma_scale;
candidates = zeros(2, numel(offsets_local)^2);
index = 0;
for dr = offsets_local
    for dtheta = offsets_local
        index = index + 1;
        candidates(:, index) = [max(cfg.channel.min_range, x_hat(1) + dr * sigma_r); ...
            wrap_angle(x_hat(2) + dtheta * sigma_theta)];
    end
end

uncertainty_points = local_points(x_hat, sigma_r, sigma_theta, cfg);
uncertainty_weights = local_weights(x_hat, uncertainty_points, sigma_r, sigma_theta);
guard_points = guard_ring_points(x_hat, sigma_r, sigma_theta, cfg);
y_average = mean(y, 2);
y_norm = max(real(y_average' * y_average), eps);
scores = -inf(1, size(candidates, 2));
beams = cell(1, size(candidates, 2));
for i = 1:size(candidates, 2)
    beams{i} = focus_beam(candidates(1, i), candidates(2, i), bits, cfg);
    main_gains = zeros(1, size(uncertainty_points, 2));
    guard_gains = zeros(1, size(guard_points, 2));
    for p = 1:size(uncertainty_points, 2)
        main_gains(p) = beam_gain(uncertainty_points(:, p), beams{i}, cfg);
    end
    for p = 1:size(guard_points, 2)
        guard_gains(p) = beam_gain(guard_points(:, p), beams{i}, cfg);
    end
    predicted_score = sum(uncertainty_weights .* main_gains) - cfg.beam.robust_penalty * ...
        (std(main_gains) + mean(guard_gains));
    measurement_score = abs(beams{i}.weights' * y_average)^2 / y_norm;
    scores(i) = (1 - cfg.beam.measurement_blend) * predicted_score + ...
        cfg.beam.measurement_blend * measurement_score;
end
[best_score, best] = max(scores);
beam = beams{best};
selection = struct('mode', 'robust_local', 'score', best_score, ...
    'candidate', candidates(:, best), 'num_candidates', numel(scores));
end

function points = local_points(x_hat, sigma_r, sigma_theta, cfg)
offsets = cfg.beam.uncertainty_offsets;
points = zeros(2, numel(offsets)^2);
index = 0;
for dr = offsets
    for dtheta = offsets
        index = index + 1;
        points(:, index) = [x_hat(1) + dr * sigma_r; x_hat(2) + dtheta * sigma_theta];
    end
end
points(1, :) = max(points(1, :), cfg.channel.min_range);
points(2, :) = wrap_angle(points(2, :));
end

function weights = local_weights(x_hat, points, sigma_r, sigma_theta)
normalized_r = (points(1, :) - x_hat(1)) ./ max(sigma_r, eps);
normalized_theta = wrap_angle(points(2, :) - x_hat(2)) ./ max(sigma_theta, eps);
weights = exp(-0.5 * (normalized_r.^2 + normalized_theta.^2));
weights = weights ./ sum(weights);
end

function points = guard_ring_points(x_hat, sigma_r, sigma_theta, cfg)
scale = cfg.beam.guard_sigma_scale;
points = [x_hat(1) + scale * sigma_r, x_hat(1) - scale * sigma_r, ...
    x_hat(1), x_hat(1), x_hat(1) + scale * sigma_r, x_hat(1) - scale * sigma_r; ...
    x_hat(2), x_hat(2), x_hat(2) + scale * sigma_theta, ...
    x_hat(2) - scale * sigma_theta, x_hat(2) + scale * sigma_theta, ...
    x_hat(2) - scale * sigma_theta];
points(1, :) = max(points(1, :), cfg.channel.min_range);
points(2, :) = wrap_angle(points(2, :));
end

function [beam, selection] = reacquire_beam(y, cfg, bits)
% Measurement-based broad codebook search after a declared loss of lock.
y_average = mean(y, 2);
y_norm = max(real(y_average' * y_average), eps);
best_metric = -inf;
beam = focus_beam(cfg.beam.reacquire_ranges(1), cfg.beam.reacquire_angles(1), bits, cfg);
for r = cfg.beam.reacquire_ranges
    for theta = cfg.beam.reacquire_angles
        candidate = focus_beam(r, theta, bits, cfg);
        metric = abs(candidate.weights' * y_average)^2 / y_norm;
        if metric > best_metric
            best_metric = metric;
            beam = candidate;
        end
    end
end
selection = struct('mode', 'reacquisition', 'score', best_metric, ...
    'candidate', beam.focus, 'num_candidates', numel(cfg.beam.reacquire_ranges) * ...
    numel(cfg.beam.reacquire_angles));
end
