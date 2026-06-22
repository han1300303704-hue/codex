function z = measurement_model(state, offsets, cfg, model_hardware)
%MEASUREMENT_MODEL Stacked real/imaginary multi-pilot array observation.

if nargin < 4
    model_hardware = true;
end
channel = spherical_channel(state(1), state(2), cfg);
base = sqrt(cfg.pilot.power) * channel.h;
y_hat = zeros(cfg.n_ant, numel(offsets));
for m = 1:numel(offsets)
    if model_hardware
        phase = state(6) + 2 * pi * state(5) * offsets(m);
        y_hat(:, m) = base * exp(1j * phase);
    else
        y_hat(:, m) = base;
    end
end
z = [real(y_hat(:)); imag(y_hat(:))];
end
