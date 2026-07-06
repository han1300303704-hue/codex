function z = measurement_model(state, offsets, cfg, model_hardware, observation_beam)
%MEASUREMENT_MODEL Stacked real/imaginary multi-pilot array observation.

if nargin < 4
    model_hardware = true;
end
if nargin < 5
    observation_beam = [];
end
channel = spherical_channel(state(1), state(2), cfg);
if isempty(observation_beam)
    effective_beam_gain = 1;
else
    effective_beam_gain = max(beam_gain(state, observation_beam, cfg), cfg.pilot.beam_gain_floor);
end
base = sqrt(cfg.pilot.power * effective_beam_gain) * channel.h;
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
