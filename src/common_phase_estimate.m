function phi = common_phase_estimate(y, predicted_state, offsets, cfg)
%COMMON_PHASE_ESTIMATE Estimate the common phase after CFO de-rotation.

a = spherical_channel(predicted_state(1), predicted_state(2), cfg).steering;
projection = zeros(1, numel(offsets));
for m = 1:numel(offsets)
    de_rotated = y(:, m) * exp(-1j * 2 * pi * predicted_state(5) * offsets(m));
    projection(m) = a' * de_rotated;
end
phi = angle(sum(projection));
end
