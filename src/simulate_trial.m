function trial = simulate_trial(cfg, scenario, seed)
%SIMULATE_TRIAL Run one closed-loop tracking realization for one scenario.

rng(seed, 'twister');
n_slots = cfg.n_slots;
truth = cfg.initial_state(:);
truth(5) = cfg.hardware.cfo_hz;
if ~scenario.physical_rf
    truth(5:6) = 0;
end

x_hat = truth + cfg.initial_error_std(:) .* randn(6, 1);
if ~scenario.model_rf
    x_hat(5:6) = 0;
end
x_hat(1) = max(x_hat(1), cfg.channel.min_range);
x_hat(2) = wrap_angle(x_hat(2));
P = cfg.initial_cov;

gain = zeros(1, n_slots);
range_error = zeros(1, n_slots);
angle_error = zeros(1, n_slots);
position_error = zeros(1, n_slots);
cfo_error = zeros(1, n_slots);
innovation = zeros(1, n_slots);
doppler_spread_hz = zeros(1, n_slots);
reacquire = false;
loss_counter = 0;
loss_events = 0;
reacquisition_count = 0;
outage_slots = 0;
beam_modes = cell(1, n_slots);
current_beam = [];

for k = 1:n_slots
    if k > 1
        truth = propagate_state(truth, cfg, true, scenario.physical_rf);
    end
    if k == 1
        n_pilots = cfg.pilot.initial_burst;
    else
        n_pilots = cfg.pilot.per_slot;
    end
    if reacquire
        observation_beam = [];
    else
        observation_beam = current_beam;
    end
    [y, offsets] = pilot_observation(truth, cfg, n_pilots, scenario.physical_rf, observation_beam);
    if k == 1
        [initial_position, ~] = coarse_position_estimate(y, x_hat, cfg);
        x_hat(1:2) = initial_position;
        P(1, 1) = cfg.initializer.range_std_m^2;
        P(2, 2) = cfg.initializer.angle_std_rad^2;
        x_pred = x_hat;
        P_pred = P;
    else
        [x_pred, P_pred] = ekf_predict(x_hat, P, cfg, scenario.model_rf);
    end

    if scenario.model_rf
        if k == 1
            if scenario.physical_rf
                x_pred(5) = coarse_cfo_estimate(y, cfg.pilot.spacing);
            else
                x_pred(5:6) = 0;
            end
            P_pred(5, 5) = min(P_pred(5, 5), cfg.initial_cov(5, 5));
        end
        if scenario.physical_rf
            phi_measured = common_phase_estimate(y, x_pred, offsets, cfg);
            x_pred(6) = phase_blend(x_pred(6), phi_measured, cfg.filter.cpe_blend);
        end
    end

    [x_hat, P, diagnostics] = ekf_update(x_pred, P_pred, y, offsets, cfg, scenario.model_rf, observation_beam);
    if scenario.model_rf && cfg.tracker.local_refinement
        [refined_position, ~] = refine_position_estimate(y, x_hat, cfg);
        x_hat(1:2) = refined_position;
        P(1, 1) = min(P(1, 1), cfg.tracker.refine_range_std_m^2);
        P(2, 2) = min(P(2, 2), cfg.tracker.refine_angle_std_rad^2);
    end
    [beam, selection] = select_beam(x_hat, P, y, offsets, cfg, ...
        scenario.quantized, scenario.robust, reacquire);
    if reacquire
        reacquisition_count = reacquisition_count + 1;
    end
    reacquire = false;
    current_beam = beam;

    gain(k) = beam_gain(truth, beam, cfg);
    range_error(k) = abs(x_hat(1) - truth(1));
    angle_error(k) = abs(wrap_angle(x_hat(2) - truth(2)));
    position_error(k) = euclidean_position_error(x_hat, truth);
    cfo_error(k) = abs(x_hat(5) - truth(5));
    innovation(k) = diagnostics.innovation_energy;
    doppler_spread_hz(k) = std(nonuniform_doppler(truth, cfg));
    beam_modes{k} = selection.mode;

    if gain(k) < cfg.beam.loss_gain_threshold
        outage_slots = outage_slots + 1;
        loss_counter = loss_counter + 1;
    else
        loss_counter = 0;
    end
    if loss_counter >= cfg.beam.loss_consecutive_slots
        loss_events = loss_events + 1;
        reacquire = true;
        loss_counter = 0;
    end
end

trial = struct();
trial.gain = gain;
trial.range_error = range_error;
trial.angle_error = angle_error;
trial.position_error = position_error;
trial.cfo_error = cfo_error;
trial.innovation = innovation;
trial.doppler_spread_hz = doppler_spread_hz;
trial.loss_events = loss_events;
trial.lost = loss_events > 0;
trial.outage_fraction = outage_slots / n_slots;
trial.reacquisition_count = reacquisition_count;
trial.beam_modes = beam_modes;
end

function phase = phase_blend(predicted, observed, blend)
phase = angle((1 - blend) * exp(1j * predicted) + blend * exp(1j * observed));
end

function error = euclidean_position_error(estimate, truth)
p_est = estimate(1) * [sin(estimate(2)); cos(estimate(2))];
p_true = truth(1) * [sin(truth(2)); cos(truth(2))];
error = norm(p_est - p_true);
end
