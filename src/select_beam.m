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
use_relaxed_quantization = use_quantization && ...
    isfield(cfg.beam, 'relaxed_quantization') && cfg.beam.relaxed_quantization && ...
    bits <= cfg.beam.relaxed_quantization_bits_max;
use_cvx_sdr = use_quantization && ...
    isfield(cfg.beam, 'cvx_sdr_quantization') && cfg.beam.cvx_sdr_quantization && ...
    bits <= cfg.beam.cvx_sdr_bits_max && exist('cvx_begin', 'file') == 2;

y_average = mean(y, 2);
y_norm = max(real(y_average' * y_average), eps);
scores = [];
beams = {};
candidate_list = [];
method_list = {};
for i = 1:size(candidates, 2)
    focus_candidate = focus_beam(candidates(1, i), candidates(2, i), bits, cfg);
    candidate_beams = {focus_candidate};
    candidate_methods = {'focus_quantized'};
    if use_relaxed_quantization
        relaxed_candidate = relaxed_quantized_beam(candidates(:, i), uncertainty_points, ...
            uncertainty_weights, bits, cfg);
        focus_center_gain = beam_gain(candidates(:, i), focus_candidate, cfg);
        relaxed_center_gain = beam_gain(candidates(:, i), relaxed_candidate, cfg);
        if relaxed_center_gain >= cfg.beam.relaxed_min_center_gain_ratio * max(focus_center_gain, eps)
            candidate_beams{end + 1} = relaxed_candidate; %#ok<AGROW>
            candidate_methods{end + 1} = 'relaxed_region'; %#ok<AGROW>
        end
    end
    for b = 1:numel(candidate_beams)
        beam_index = numel(beams) + 1;
        beams{beam_index} = candidate_beams{b}; %#ok<AGROW>
        candidate_list(:, beam_index) = candidates(:, i); %#ok<AGROW>
        method_list{beam_index} = candidate_methods{b}; %#ok<AGROW>
        scores(beam_index) = score_beam(candidate_beams{b}, uncertainty_points, uncertainty_weights, ...
            guard_points, y_average, y_norm, cfg); %#ok<AGROW>
    end
end
[best_score, best] = max(scores);
beam = beams{best};
best_candidate = candidate_list(:, best);
best_method = method_list{best};
if use_cvx_sdr
    cvx_candidate = cvx_sdr_quantized_beam(best_candidate, uncertainty_points, uncertainty_weights, bits, cfg);
    reference_gain = beam_gain(best_candidate, beam, cfg);
    cvx_center_gain = beam_gain(best_candidate, cvx_candidate, cfg);
    if cvx_center_gain >= cfg.beam.relaxed_min_center_gain_ratio * max(reference_gain, eps)
        cvx_score = score_beam(cvx_candidate, uncertainty_points, uncertainty_weights, ...
            guard_points, y_average, y_norm, cfg);
        scores(end + 1) = cvx_score; %#ok<AGROW>
        if cvx_score > best_score
            best_score = cvx_score;
            beam = cvx_candidate;
            best_method = cvx_candidate.quantization.method;
        end
    end
end
selection = struct('mode', ['robust_local_' best_method], 'score', best_score, ...
    'candidate', best_candidate, 'num_candidates', numel(scores));
end

function score = score_beam(beam, uncertainty_points, uncertainty_weights, guard_points, y_average, y_norm, cfg)
    main_gains = zeros(1, size(uncertainty_points, 2));
    guard_gains = zeros(1, size(guard_points, 2));
    for p = 1:size(uncertainty_points, 2)
        main_gains(p) = beam_gain(uncertainty_points(:, p), beam, cfg);
    end
    for p = 1:size(guard_points, 2)
        guard_gains(p) = beam_gain(guard_points(:, p), beam, cfg);
    end
    predicted_score = sum(uncertainty_weights .* main_gains) - cfg.beam.robust_penalty * ...
        (std(main_gains) + mean(guard_gains));
    measurement_score = abs(beam.weights' * y_average)^2 / y_norm;
    score = (1 - cfg.beam.measurement_blend) * predicted_score + ...
        cfg.beam.measurement_blend * measurement_score;
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
