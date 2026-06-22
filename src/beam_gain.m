function gain = beam_gain(state_or_position, beam, cfg)
%BEAM_GAIN Normalized receive power at a state/position for a given beam.

channel = spherical_channel(state_or_position(1), state_or_position(2), cfg);
h = channel.h;
gain = abs(h' * beam.weights)^2 / max(real(h' * h), eps);
gain = min(max(real(gain), 0), 1.05);
end
