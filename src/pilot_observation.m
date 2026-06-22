function [y, offsets] = pilot_observation(true_state, cfg, n_pilots, include_hardware)
%PILOT_OBSERVATION Generate N-by-M uplink pilot observations in one slot.
% Phase noise is common to the array but evolves within the pilot mini-burst.

if nargin < 3 || isempty(n_pilots)
    n_pilots = cfg.pilot.per_slot;
end
if nargin < 4
    include_hardware = true;
end

offsets = (0:n_pilots-1) * cfg.pilot.spacing;
channel = spherical_channel(true_state(1), true_state(2), cfg);
signal = sqrt(cfg.pilot.power) * channel.h;
noise_variance = cfg.pilot.power / (10^(cfg.tx_snr_db / 10));
y = zeros(cfg.n_ant, n_pilots);

phase = true_state(6);
for m = 1:n_pilots
    if include_hardware
        if m > 1
            dt = offsets(m) - offsets(m-1);
            phase = wrap_angle(phase + 2 * pi * true_state(5) * dt + ...
                sqrt(2 * pi * cfg.hardware.phase_noise_linewidth_hz * dt) * randn);
        end
        impairment = exp(1j * phase);
    else
        impairment = 1;
    end
    noise = sqrt(noise_variance / 2) * (randn(cfg.n_ant, 1) + 1j * randn(cfg.n_ant, 1));
    y(:, m) = signal * impairment + noise;
end
end
