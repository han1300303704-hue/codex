function beam = focus_beam(r, theta, bits, cfg)
%FOCUS_BEAM Spherical focusing beam with optional B-bit phase quantization.

channel = spherical_channel(r, theta, cfg);
continuous = channel.steering;
if isinf(bits)
    weights = continuous;
    quantized_phase = angle(weights);
    optimization = struct('method', 'continuous', 'offset_rad', 0, 'gain_metric', NaN);
else
    [weights, quantized_phase, optimization] = quantized_focus_weights(continuous, channel.h, bits, cfg);
end
beam = struct('weights', weights, 'continuous_phase', angle(continuous), ...
    'quantized_phase', quantized_phase, 'bits', bits, 'focus', [r; theta], ...
    'quantization', optimization);
end

function [weights, quantized_phase, optimization] = quantized_focus_weights(continuous, h, bits, cfg)
levels = 2^bits;
step = 2 * pi / levels;
continuous_phase = mod(angle(continuous), 2 * pi);

best_phase = nearest_phase(continuous_phase, step);
best_metric = phase_metric(best_phase, h, cfg.n_ant);
best_offset = 0;
method = 'nearest';

if isfield(cfg.beam, 'quantization_offset_search') && cfg.beam.quantization_offset_search && ...
        bits <= cfg.beam.quantization_optimize_bits_max
    n_trials = max(1, round(cfg.beam.quantization_offset_trials));
    offsets = (0:n_trials-1) * step / n_trials;
    for i = 1:numel(offsets)
        candidate_phase = nearest_phase(continuous_phase + offsets(i), step);
        candidate_metric = phase_metric(candidate_phase, h, cfg.n_ant);
        if candidate_metric > best_metric
            best_metric = candidate_metric;
            best_phase = candidate_phase;
            best_offset = offsets(i);
            method = 'offset_search';
        end
    end
end

greedy_passes = 0;
if isfield(cfg.beam, 'quantization_greedy_passes')
    greedy_passes = max(0, round(cfg.beam.quantization_greedy_passes));
end
if greedy_passes > 0 && bits <= cfg.beam.quantization_optimize_bits_max
    [best_phase, best_metric] = greedy_neighbor_refinement(best_phase, h, step, levels, cfg.n_ant, greedy_passes);
    method = [method '+greedy'];
end

quantized_phase = mod(best_phase, 2 * pi);
weights = exp(1j * quantized_phase) ./ sqrt(cfg.n_ant);
optimization = struct('method', method, 'offset_rad', best_offset, 'gain_metric', best_metric);
end

function phase = nearest_phase(phase_in, step)
phase = mod(step * round(mod(phase_in, 2 * pi) / step), 2 * pi);
end

function metric = phase_metric(phase, h, n_ant)
weights = exp(1j * phase) ./ sqrt(n_ant);
metric = abs(h' * weights)^2;
end

function [phase, metric] = greedy_neighbor_refinement(phase, h, step, levels, n_ant, passes)
metric = phase_metric(phase, h, n_ant);
for pass = 1:passes
    improved = false;
    order = 1:numel(phase);
    for n = order
        current = phase(n);
        best_local_phase = current;
        best_local_metric = metric;
        for delta = [-1 1]
            trial_phase = phase;
            trial_phase(n) = mod(current + delta * step, 2 * pi);
            trial_metric = phase_metric(trial_phase, h, n_ant);
            if trial_metric > best_local_metric
                best_local_metric = trial_metric;
                best_local_phase = trial_phase(n);
            end
        end
        if best_local_metric > metric
            phase(n) = best_local_phase;
            metric = best_local_metric;
            improved = true;
        end
    end
    if ~improved || levels <= 1
        break;
    end
end
end
