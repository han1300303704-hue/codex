function result = run_monte_carlo(cfg, scenario)
%RUN_MONTE_CARLO Aggregate time-domain tracking statistics over trials.

n_mc = cfg.n_mc;
n_slots = cfg.n_slots;
gain_sum = zeros(1, n_slots);
range_sq_sum = zeros(1, n_slots);
angle_sq_sum = zeros(1, n_slots);
position_sq_sum = zeros(1, n_slots);
cfo_sq_sum = zeros(1, n_slots);
doppler_spread_sum = zeros(1, n_slots);
lost_count = 0;
loss_event_count = 0;
reacquisition_count = 0;

for mc = 1:n_mc
    % Common random numbers keep ablation comparisons statistically fair.
    trial_seed = cfg.seed + mc;
    trial = simulate_trial(cfg, scenario, trial_seed);
    gain_sum = gain_sum + trial.gain;
    range_sq_sum = range_sq_sum + trial.range_error.^2;
    angle_sq_sum = angle_sq_sum + trial.angle_error.^2;
    position_sq_sum = position_sq_sum + trial.position_error.^2;
    cfo_sq_sum = cfo_sq_sum + trial.cfo_error.^2;
    doppler_spread_sum = doppler_spread_sum + trial.doppler_spread_hz;
    lost_count = lost_count + double(trial.lost);
    loss_event_count = loss_event_count + trial.loss_events;
    reacquisition_count = reacquisition_count + trial.reacquisition_count;
end

result = struct();
result.name = scenario.name;
result.label = scenario.label;
result.time_s = (0:n_slots-1) * cfg.track_period;
result.mean_gain = gain_sum / n_mc;
result.range_rmse_m = sqrt(range_sq_sum / n_mc);
result.angle_rmse_rad = sqrt(angle_sq_sum / n_mc);
result.position_rmse_m = sqrt(position_sq_sum / n_mc);
result.cfo_rmse_hz = sqrt(cfo_sq_sum / n_mc);
result.mean_doppler_spread_hz = doppler_spread_sum / n_mc;
result.loss_probability = lost_count / n_mc;
result.mean_loss_events = loss_event_count / n_mc;
result.mean_reacquisitions = reacquisition_count / n_mc;
result.n_mc = n_mc;
end
