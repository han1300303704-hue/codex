function beam = cvx_sdr_quantized_beam(focus, points, weights, bits, cfg)
%CVX_SDR_QUANTIZED_BEAM CVX semidefinite relaxation for finite-phase beams.
%   For small arrays this solves the element-level SDR
%       maximize trace(R Q), diag(Q)=1, Q >= 0
%   and projects randomized relaxed vectors to the B-bit phase alphabet.
%
%   For large arrays, e.g. N=1024, the full SDP is not practical online.
%   The function then optimizes block-wise residual phase corrections on top
%   of the ordinary quantized focus beam.  This still uses a genuine CVX SDP,
%   but with dimension equal to the configured block count.

if isinf(bits)
    beam = focus_beam(focus(1), focus(2), inf, cfg);
    return;
end
if exist('cvx_begin', 'file') ~= 2
    beam = relaxed_quantized_beam(focus, points, weights, bits, cfg);
    beam.quantization.method = 'cvx_unavailable_fallback';
    return;
end

weights = weights(:).';
weights = weights ./ max(sum(weights), eps);
responses = normalized_responses(points, cfg);

if cfg.n_ant <= cfg.beam.cvx_sdr_max_full_ant
    beam = full_element_sdr(focus, responses, weights, bits, cfg);
else
    beam = block_residual_sdr(focus, responses, weights, bits, cfg);
end
end

function beam = full_element_sdr(focus, responses, weights, bits, cfg)
n = cfg.n_ant;
R_points = zeros(n, n, numel(weights));
for i = 1:numel(weights)
    a = responses(:, i);
    R_points(:, :, i) = a * a';
end
Q = solve_sdr(R_points, weights, cfg);

candidates = randomized_phases(Q, bits, cfg);
focus_beam_candidate = focus_beam(focus(1), focus(2), bits, cfg);
candidates(:, end + 1) = focus_beam_candidate.quantized_phase;

[phase, metric] = best_phase_candidate(candidates, responses, weights, cfg.n_ant, cfg);
beam = make_quantized_beam(phase, focus, bits, metric, ['cvx_sdr_full_element_' cfg.beam.cvx_sdr_objective], cfg);
end

function beam = block_residual_sdr(focus, responses, weights, bits, cfg)
base = focus_beam(focus(1), focus(2), bits, cfg);
groups = antenna_groups(cfg.n_ant, cfg.beam.cvx_sdr_block_count);
g_count = max(groups);

D = zeros(g_count, size(responses, 2));
for i = 1:size(responses, 2)
    a = responses(:, i);
    for g = 1:g_count
        idx = groups == g;
        D(g, i) = sum(conj(a(idx)) .* base.weights(idx));
    end
end

R_points = zeros(g_count, g_count, numel(weights));
for i = 1:numel(weights)
    d = D(:, i);
    R_points(:, :, i) = conj(d) * d.';
end
Q = solve_sdr(R_points, weights, cfg);

correction_candidates = randomized_phases(Q, bits, cfg);
correction_candidates(:, end + 1) = zeros(g_count, 1);

best_metric = -inf;
best_phase = base.quantized_phase;
for i = 1:size(correction_candidates, 2)
    phase = mod(base.quantized_phase + correction_candidates(groups, i), 2 * pi);
    metric = objective_metric(phase, responses, weights, cfg.n_ant, cfg);
    if metric > best_metric
        best_metric = metric;
        best_phase = phase;
    end
end
beam = make_quantized_beam(best_phase, focus, bits, best_metric, ['cvx_sdr_block_residual_' cfg.beam.cvx_sdr_objective], cfg);
end

function Q = solve_sdr(R_points, weights, cfg)
n = size(R_points, 1);
n_points = size(R_points, 3);
if isfield(cfg.beam, 'cvx_sdr_solver') && ~isempty(cfg.beam.cvx_sdr_solver)
    cvx_solver(cfg.beam.cvx_sdr_solver);
end
objective = 'average';
if isfield(cfg.beam, 'cvx_sdr_objective')
    objective = lower(cfg.beam.cvx_sdr_objective);
end
cvx_begin sdp quiet
    variable Q(n, n) hermitian semidefinite
    if strcmp(objective, 'maxmin')
        variable t
        maximize(t)
        subject to
            diag(Q) == ones(n, 1);
            for p = 1:n_points
                real(trace(R_points(:, :, p) * Q)) >= t;
            end
    else
        R = zeros(n, n);
        for p = 1:n_points
            R = R + weights(p) * R_points(:, :, p);
        end
        R = (R + R') / 2;
        maximize(real(trace(R * Q)))
        subject to
            diag(Q) == ones(n, 1);
    end
cvx_end
Q = full((Q + Q') / 2);
if ~strcmpi(cvx_status, 'Solved') && ~strcmpi(cvx_status, 'Inaccurate/Solved')
    warning('cvx_sdr_quantized_beam:CVXStatus', 'CVX returned status: %s', cvx_status);
end
end

function candidates = randomized_phases(Q, bits, cfg)
n = size(Q, 1);
levels = 2^bits;
step = 2 * pi / levels;
n_random = max(0, round(cfg.beam.cvx_sdr_randomizations));
n_offsets = 1;
if isfield(cfg.beam, 'quantization_offset_trials')
    n_offsets = max(1, round(cfg.beam.quantization_offset_trials));
end
offsets = (0:n_offsets-1) * step / n_offsets;

[V, D] = eig((Q + Q') / 2);
eigvals = max(real(diag(D)), 0);
rootQ = V * diag(sqrt(eigvals));

candidates = zeros(n, (n_random + 1) * n_offsets + 1);
col = 0;
for k = 1:n_random
    z = rootQ * (randn(n, 1) + 1j * randn(n, 1)) / sqrt(2);
    [candidates, col] = add_offset_candidates(candidates, col, angle(z), offsets, step);
end
[~, dominant] = max(eigvals);
[candidates, col] = add_offset_candidates(candidates, col, angle(V(:, dominant)), offsets, step);
candidates(:, col + 1) = zeros(n, 1);
end

function [candidates, col] = add_offset_candidates(candidates, col, phase, offsets, step)
for i = 1:numel(offsets)
    col = col + 1;
    candidates(:, col) = quantize_phase(phase + offsets(i), step);
end
end

function [best_phase, best_metric] = best_phase_candidate(candidates, responses, weights, n_ant, cfg)
best_metric = -inf;
best_phase = candidates(:, 1);
for i = 1:size(candidates, 2)
    metric = objective_metric(candidates(:, i), responses, weights, n_ant, cfg);
    if metric > best_metric
        best_metric = metric;
        best_phase = candidates(:, i);
    end
end
end

function beam = make_quantized_beam(phase, focus, bits, metric, method, cfg)
phase = mod(phase(:), 2 * pi);
beam = struct('weights', exp(1j * phase) ./ sqrt(cfg.n_ant), ...
    'continuous_phase', phase, ...
    'quantized_phase', phase, ...
    'bits', bits, ...
    'focus', focus(:), ...
    'quantization', struct('method', method, 'offset_rad', NaN, 'gain_metric', metric));
end

function responses = normalized_responses(points, cfg)
responses = zeros(cfg.n_ant, size(points, 2));
for i = 1:size(points, 2)
    channel = spherical_channel(points(1, i), points(2, i), cfg);
    responses(:, i) = channel.h ./ max(norm(channel.h), eps);
end
end

function groups = antenna_groups(n_ant, requested_groups)
g_count = min(n_ant, max(1, round(requested_groups)));
groups = ceil((1:n_ant).' * g_count / n_ant);
groups = min(max(groups, 1), g_count);
end

function phase = quantize_phase(phase, step)
phase = mod(step * round(mod(phase, 2 * pi) / step), 2 * pi);
end

function metric = region_metric(phase, responses, weights, n_ant)
w = exp(1j * phase(:)) ./ sqrt(n_ant);
metric = 0;
for i = 1:numel(weights)
    metric = metric + weights(i) * abs(responses(:, i)' * w)^2;
end
end

function metric = objective_metric(phase, responses, weights, n_ant, cfg)
objective = 'average';
if isfield(cfg.beam, 'cvx_sdr_objective')
    objective = lower(cfg.beam.cvx_sdr_objective);
end
w = exp(1j * phase(:)) ./ sqrt(n_ant);
point_gains = zeros(1, numel(weights));
for i = 1:numel(weights)
    point_gains(i) = abs(responses(:, i)' * w)^2;
end
if strcmp(objective, 'maxmin')
    metric = min(point_gains);
else
    metric = sum(weights .* point_gains);
end
end
