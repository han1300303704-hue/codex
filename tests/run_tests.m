function run_tests()
%RUN_TESTS Deterministic smoke tests for the MATLAB simulation framework.

root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'src'));
fprintf('Running near-field beam-tracking tests...\n');
test_focus_peak();
test_quantizer_levels();
test_ideal_tracking_convergence();
test_joint_improves_over_uncompensated();
test_result_artifacts();
fprintf('All tests passed.\n');
end

function cfg = small_cfg()
cfg = default_config(struct( ...
    'n_ant', 32, 'n_slots', 14, 'n_mc', 1, 'figure_visible', 'off', ...
    'initial_state', [2.0; deg2rad(15); 0.5; 0.4; 800; 0], ...
    'initial_error_std', [0.12; deg2rad(0.6); 0.08; 0.08; 80; deg2rad(8)], ...
    'scan', struct('phase_bits', 2, 'phase_linewidth_hz', 30, 'n_mc', 1)));
cfg.hardware.cfo_hz = 800;
cfg.hardware.phase_noise_linewidth_hz = 30;
cfg.pilot.initial_burst = 8;
end

function test_focus_peak()
cfg = small_cfg();
state = cfg.initial_state;
beam = focus_beam(state(1), state(2), inf, cfg);
gain_at_focus = beam_gain(state, beam, cfg);
offset_state = state;
offset_state(2) = offset_state(2) + deg2rad(8);
gain_off_focus = beam_gain(offset_state, beam, cfg);
assert(gain_at_focus > 0.98, 'Continuous focus beam must peak at its target.');
assert(gain_at_focus > gain_off_focus, 'Focus beam must reject an off-target angle.');
end

function test_quantizer_levels()
cfg = small_cfg();
bits = 3;
beam = focus_beam(cfg.initial_state(1), cfg.initial_state(2), bits, cfg);
step = 2 * pi / 2^bits;
phase_index = beam.quantized_phase / step;
assert(max(abs(phase_index - round(phase_index))) < 1e-10, ...
    'Quantized beam phases must lie on valid B-bit levels.');
assert(abs(norm(beam.weights) - 1) < 1e-10, 'Beam weights must have unit norm.');
end

function test_ideal_tracking_convergence()
cfg = small_cfg();
scenarios = make_scenarios();
trial = simulate_trial(cfg, scenarios(1), cfg.seed + 1);
assert(all(isfinite(trial.position_error)), 'Ideal trial must produce finite tracking errors.');
assert(mean(trial.gain(end-2:end)) > 0.45, ...
    'Ideal tracking must retain a usable focus gain after convergence.');
end

function test_joint_improves_over_uncompensated()
cfg = small_cfg();
cfg.hardware.cfo_hz = 2500;
cfg.hardware.phase_noise_linewidth_hz = 300;
scenarios = make_scenarios();
uncompensated = simulate_trial(cfg, scenarios(2), cfg.seed + 7);
joint = simulate_trial(cfg, scenarios(5), cfg.seed + 7);
assert(mean(joint.gain(end-4:end)) > mean(uncompensated.gain(end-4:end)), ...
    'Joint compensation should improve fixed-seed late-slot gain over no compensation.');
end

function test_result_artifacts()
cfg = small_cfg();
cfg.n_slots = 5;
cfg.output_dir = tempname;
mkdir(cfg.output_dir);
run_suite(cfg);
required = {'tracking_results.mat', 'beam_gain.fig', 'beam_gain.png', ...
    'tracking_error.fig', 'tracking_error.png', 'lock_loss_probability.fig', ...
    'lock_loss_probability.png', 'robustness_scans.fig', 'robustness_scans.png'};
for i = 1:numel(required)
    assert(exist(fullfile(cfg.output_dir, required{i}), 'file') == 2, ...
        ['Missing expected result artifact: ' required{i}]);
end
end
