function run_tests()
%RUN_TESTS Deterministic smoke tests for the MATLAB simulation framework.

root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'src'));
fprintf('Running near-field beam-tracking tests...\n');
test_focus_peak();
test_quantizer_levels();
test_offset_quantizer_not_worse();
test_relaxed_quantizer_levels();
test_cvx_sdr_quantizer_if_available();
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

function test_offset_quantizer_not_worse()
cfg = small_cfg();
state = cfg.initial_state;
nearest_cfg = cfg;
nearest_cfg.beam.quantization_offset_search = false;
optimized_cfg = cfg;
optimized_cfg.beam.quantization_offset_search = true;
for bits = 1:2
    nearest = focus_beam(state(1), state(2), bits, nearest_cfg);
    optimized = focus_beam(state(1), state(2), bits, optimized_cfg);
    assert(beam_gain(state, optimized, optimized_cfg) + 1e-12 >= ...
        beam_gain(state, nearest, nearest_cfg), ...
        'Offset-search quantization must not reduce focus gain.');
end
end

function test_relaxed_quantizer_levels()
cfg = small_cfg();
state = cfg.initial_state;
points = [state(1), state(1) + 0.05, state(1) - 0.05; ...
    state(2), state(2) + deg2rad(0.2), state(2) - deg2rad(0.2)];
weights = [0.5, 0.25, 0.25];
bits = 2;
beam = relaxed_quantized_beam(state(1:2), points, weights, bits, cfg);
step = 2 * pi / 2^bits;
phase_index = beam.quantized_phase / step;
assert(max(abs(phase_index - round(phase_index))) < 1e-10, ...
    'Relaxed quantized beam phases must lie on valid B-bit levels.');
assert(abs(norm(beam.weights) - 1) < 1e-10, 'Relaxed beam weights must have unit norm.');
end

function test_cvx_sdr_quantizer_if_available()
if exist('cvx_begin', 'file') ~= 2
    fprintf('Skipping CVX-SDR quantizer test because CVX is not on the MATLAB path.\n');
    return;
end
cfg = default_config(struct('n_ant', 8, 'n_slots', 4, 'n_mc', 1, 'figure_visible', 'off'));
cfg.beam.cvx_sdr_max_full_ant = 16;
cfg.beam.cvx_sdr_randomizations = 4;
state = [2.0; deg2rad(15); 0; 0; 0; 0];
points = [state(1), state(1) + 0.03, state(1) - 0.03; ...
    state(2), state(2) + deg2rad(0.1), state(2) - deg2rad(0.1)];
weights = [0.5, 0.25, 0.25];
bits = 2;
beam = cvx_sdr_quantized_beam(state(1:2), points, weights, bits, cfg);
step = 2 * pi / 2^bits;
phase_index = beam.quantized_phase / step;
assert(max(abs(phase_index - round(phase_index))) < 1e-10, ...
    'CVX-SDR quantized beam phases must lie on valid B-bit levels.');
assert(abs(norm(beam.weights) - 1) < 1e-10, 'CVX-SDR beam weights must have unit norm.');
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
joint = simulate_trial(cfg, get_scenario(scenarios, 'joint'), cfg.seed + 7);
assert(mean(joint.position_error(end-4:end)) < mean(uncompensated.position_error(end-4:end)), ...
    'Joint compensation should improve fixed-seed late-slot position tracking over no compensation.');
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

function scenario = get_scenario(scenarios, name)
idx = strcmp({scenarios.name}, name);
assert(any(idx), ['Missing scenario: ' name]);
scenario = scenarios(find(idx, 1, 'first'));
end
