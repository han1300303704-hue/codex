function beam = focus_beam(r, theta, bits, cfg)
%FOCUS_BEAM Spherical focusing beam with optional B-bit phase quantization.

channel = spherical_channel(r, theta, cfg);
continuous = channel.steering;
if isinf(bits)
    weights = continuous;
    quantized_phase = angle(weights);
else
    levels = 2^bits;
    step = 2 * pi / levels;
    continuous_phase = mod(angle(continuous), 2 * pi);
    quantized_phase = mod(step * round(continuous_phase / step), 2 * pi);
    weights = exp(1j * quantized_phase) ./ sqrt(cfg.n_ant);
end
beam = struct('weights', weights, 'continuous_phase', angle(continuous), ...
    'quantized_phase', quantized_phase, 'bits', bits, 'focus', [r; theta]);
end
