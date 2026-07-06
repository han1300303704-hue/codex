function [x_post, P_post, diagnostics] = ekf_update(x_pred, P_pred, y, offsets, cfg, model_hardware, observation_beam)
%EKF_UPDATE Information-form EKF correction for high-dimensional pilots.
% Using the information form avoids inversion of an (2NM)-by-(2NM) matrix.

if nargin < 7
    observation_beam = [];
end
z = [real(y(:)); imag(y(:))];
z_hat = measurement_model(x_pred, offsets, cfg, model_hardware, observation_beam);
residual = z - z_hat;
n_state = numel(x_pred);
H = zeros(numel(z), n_state);
steps = cfg.filter.jacobian_step(:);

for d = 1:n_state
    perturb = zeros(n_state, 1);
    perturb(d) = steps(d);
    h_plus = measurement_model(x_pred + perturb, offsets, cfg, model_hardware, observation_beam);
    h_minus = measurement_model(x_pred - perturb, offsets, cfg, model_hardware, observation_beam);
    H(:, d) = (h_plus - h_minus) / (2 * steps(d));
end

if ~model_hardware
    H(:, 5:6) = 0;
end

noise_variance = cfg.pilot.power / (10^(cfg.tx_snr_db / 10));
inv_r = 2 / noise_variance;
prior_information = pinv(P_pred);
information = prior_information + inv_r * (H' * H);
P_candidate = stabilize_covariance(pinv(information), cfg.filter.min_cov_eig);
innovation_energy = real(inv_r * (residual' * residual));
increment = P_candidate * (inv_r * (H' * residual));

% Limit catastrophic corrections while preserving normal EKF behaviour.
if innovation_energy > cfg.filter.innovation_gate * numel(residual)
    increment = 0.35 * increment;
end
x_post = x_pred + increment;
x_post(1) = max(x_post(1), cfg.channel.min_range);
x_post(2) = wrap_angle(x_post(2));
x_post(6) = wrap_angle(x_post(6));
if ~model_hardware
    x_post(5:6) = 0;
end
P_post = P_candidate;
diagnostics = struct('innovation_energy', innovation_energy, ...
    'residual_norm', norm(residual), 'jacobian', H);
end
