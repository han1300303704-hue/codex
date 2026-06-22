function [x_pred, P_pred] = ekf_predict(x, P, cfg, model_hardware)
%EKF_PREDICT Nonlinear state prediction with numerical transition Jacobian.

if nargin < 4
    model_hardware = true;
end
x_pred = propagate_state(x, cfg, false, model_hardware);
n_state = numel(x);
F = zeros(n_state);
steps = cfg.filter.jacobian_step(:);
for d = 1:n_state
    perturb = zeros(n_state, 1);
    perturb(d) = steps(d);
    plus = propagate_state(x + perturb, cfg, false, model_hardware);
    minus = propagate_state(x - perturb, cfg, false, model_hardware);
    difference = plus - minus;
    difference(2) = wrap_angle(difference(2));
    difference(6) = wrap_angle(difference(6));
    F(:, d) = difference / (2 * steps(d));
end

q = cfg.filter.process_std(:).^2;
if ~model_hardware
    q(5:6) = cfg.filter.min_cov_eig;
end
P_pred = F * P * F' + diag(q);
P_pred = stabilize_covariance(P_pred, cfg.filter.min_cov_eig);
end
