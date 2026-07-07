function beam = relaxed_quantized_beam(focus, points, weights, bits, cfg)
%RELAXED_QUANTIZED_BEAM SDR-style robust finite-phase beam.
%   The finite-bit problem
%       maximize sum_i p_i |h_i^H w|^2,  w_n in finite phase alphabet
%   is first relaxed to a continuous unit-modulus/eigenvector surrogate.
%   The relaxed vector is then projected back to the B-bit phase alphabet
%   with a small global-offset search.  This keeps the implementation pure
%   MATLAB and toolbox-free while capturing the main SDR/codebook-design
%   idea: optimize a region objective instead of quantizing a single focus.

if isinf(bits)
    beam = focus_beam(focus(1), focus(2), inf, cfg);
    return;
end

weights = weights(:).';
weights = weights ./ max(sum(weights), eps);
responses = normalized_responses(points, cfg);

v = spherical_channel(focus(1), focus(2), cfg).steering;
v = v ./ max(norm(v), eps);
iterations = max(1, round(cfg.beam.relaxed_power_iterations));
for iter = 1:iterations
    v_next = apply_region_covariance(v, responses, weights);
    if norm(v_next) < eps
        break;
    end
    v = v_next ./ norm(v_next);
end

[phase, metric, offset] = quantize_relaxed_phase(angle(v), responses, weights, bits, cfg);
beam = struct('weights', exp(1j * phase) ./ sqrt(cfg.n_ant), ...
    'continuous_phase', angle(v), ...
    'quantized_phase', phase, ...
    'bits', bits, ...
    'focus', focus(:), ...
    'quantization', struct('method', 'relaxed_region_projection', ...
    'offset_rad', offset, 'gain_metric', metric));
end

function responses = normalized_responses(points, cfg)
responses = zeros(cfg.n_ant, size(points, 2));
for i = 1:size(points, 2)
    channel = spherical_channel(points(1, i), points(2, i), cfg);
    responses(:, i) = channel.h ./ max(norm(channel.h), eps);
end
end

function y = apply_region_covariance(v, responses, weights)
y = zeros(size(v));
for i = 1:numel(weights)
    a = responses(:, i);
    y = y + weights(i) * a * (a' * v);
end
end

function [best_phase, best_metric, best_offset] = quantize_relaxed_phase(phase0, responses, weights, bits, cfg)
levels = 2^bits;
step = 2 * pi / levels;
phase0 = mod(phase0(:), 2 * pi);

n_trials = 1;
if isfield(cfg.beam, 'quantization_offset_trials')
    n_trials = max(1, round(cfg.beam.quantization_offset_trials));
end
offsets = (0:n_trials-1) * step / n_trials;

best_metric = -inf;
best_phase = zeros(size(phase0));
best_offset = 0;
for i = 1:numel(offsets)
    candidate = mod(step * round(mod(phase0 + offsets(i), 2 * pi) / step), 2 * pi);
    metric = region_metric(candidate, responses, weights, cfg.n_ant);
    if metric > best_metric
        best_metric = metric;
        best_phase = candidate;
        best_offset = offsets(i);
    end
end
end

function metric = region_metric(phase, responses, weights, n_ant)
w = exp(1j * phase) ./ sqrt(n_ant);
metric = 0;
for i = 1:numel(weights)
    metric = metric + weights(i) * abs(responses(:, i)' * w)^2;
end
end
