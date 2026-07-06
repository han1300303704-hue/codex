function cfg = default_config(overrides)
%DEFAULT_CONFIG Reproducible parameters for near-field beam tracking.
%   CFG = DEFAULT_CONFIG() returns a pure-MATLAB configuration.  Pass a
%   structure to override top-level fields, e.g. DEFAULT_CONFIG(struct('n_mc',5)).

if nargin < 1
    overrides = struct();
end

cfg = struct();
cfg.seed = 20260623;
cfg.c = 299792458;
cfg.fc = 300e9;
cfg.lambda = cfg.c / cfg.fc;
cfg.n_ant = 1024;
cfg.spacing = cfg.lambda / 2;
cfg.tx_snr_db = 10;
cfg.track_period = 50e-6;
cfg.n_slots = 200;
cfg.n_mc = 50;
cfg.figure_visible = 'off';
cfg.output_dir = fullfile(pwd, 'results');

cfg.initial_state = [1.5; deg2rad(20); 2; 3; 2e3; 0];
cfg.initial_error_std = [0.20; deg2rad(1.5); 0.35; 0.35; 350; deg2rad(20)];
cfg.initial_cov = diag(cfg.initial_error_std .^ 2);

cfg.channel = struct();
cfg.channel.path_gain_at_1m = 1;
cfg.channel.min_range = 1.0;

cfg.motion = struct();
cfg.motion.accel_std = 0.08;              % m/s^2 per tracking slot
cfg.motion.cfo_rw_std = 3.0;              % Hz per tracking slot

cfg.hardware = struct();
cfg.hardware.cfo_hz = 2e3;
cfg.hardware.phase_noise_linewidth_hz = 100;
cfg.hardware.phase_shifter_bits = 3;

cfg.pilot = struct();
cfg.pilot.per_slot = 4;
cfg.pilot.initial_burst = 32;
cfg.pilot.spacing = cfg.track_period / 8;
cfg.pilot.power = 1;
cfg.pilot.beam_gain_floor = 1e-3;

cfg.initializer = struct();
cfg.initializer.range_half_width_m = 0.8;
cfg.initializer.range_step_m = 0.01;
cfg.initializer.angle_half_width_rad = deg2rad(6);
cfg.initializer.angle_step_rad = deg2rad(0.05);
cfg.initializer.range_std_m = 0.25;
cfg.initializer.angle_std_rad = deg2rad(0.5);

cfg.tracker = struct();
cfg.tracker.local_refinement = true;
cfg.tracker.refine_range_half_width_m = 0.12;
cfg.tracker.refine_range_step_m = 0.005;
cfg.tracker.refine_angle_half_width_rad = deg2rad(0.35);
cfg.tracker.refine_angle_step_rad = deg2rad(0.01);
cfg.tracker.refine_range_std_m = 0.03;
cfg.tracker.refine_angle_std_rad = deg2rad(0.04);

cfg.filter = struct();
cfg.filter.process_std = [0.015; deg2rad(0.02); 0.10; 0.10; 8; deg2rad(3)];
cfg.filter.jacobian_step = [1e-3; 1e-5; 1e-3; 1e-3; 0.5; 1e-4];
cfg.filter.min_cov_eig = 1e-10;
cfg.filter.cpe_blend = 0.30;
cfg.filter.innovation_gate = 60;

cfg.beam = struct();
cfg.beam.local_sigma_scale = 1.5;
cfg.beam.local_offsets = -2:2;
cfg.beam.uncertainty_offsets = [-1 0 1];
cfg.beam.robust_min_range_std_m = 0.60;
cfg.beam.robust_min_angle_std_rad = deg2rad(1.0);
cfg.beam.guard_sigma_scale = 4.0;
cfg.beam.robust_penalty = 0.25;
cfg.beam.measurement_blend = 0.0;
cfg.beam.loss_gain_threshold = 0.5;
cfg.beam.loss_consecutive_slots = 3;
cfg.beam.reacquire_ranges = linspace(0.8, 4.0, 9);
cfg.beam.reacquire_angles = deg2rad(-45:3:45);

cfg.scan = struct();
cfg.scan.phase_bits = 1:4;
cfg.scan.phase_linewidth_hz = [10 30 100 300 1000];
cfg.scan.n_mc = min(30, cfg.n_mc);

cfg.stress = struct();
cfg.stress.enabled = true;
cfg.stress.n_mc = min(20, cfg.n_mc);
cfg.stress.phase_shifter_bits = 2;
cfg.stress.tx_snr_db = 0;
cfg.stress.phase_noise_linewidth_hz = 300;
cfg.stress.pilot_per_slot = 1;
cfg.stress.initial_state = [10; deg2rad(20); 2; 12; 2e3; 0];
cfg.stress.initial_error_std = [1.00; deg2rad(4.0); 0.80; 1.00; 500; deg2rad(30)];

cfg = merge_struct(cfg, overrides);
cfg.lambda = cfg.c / cfg.fc;
cfg.spacing = cfg.lambda / 2;
if isfield(overrides, 'n_mc') && (~isfield(overrides, 'scan') || ~isfield(overrides.scan, 'n_mc'))
    cfg.scan.n_mc = min(30, cfg.n_mc);
    cfg.stress.n_mc = min(20, cfg.n_mc);
end
if isfield(overrides, 'initial_error_std') && ~isfield(overrides, 'initial_cov')
    cfg.initial_cov = diag(cfg.initial_error_std(:) .^ 2);
end
if ~isfield(overrides, 'pilot') || ~isfield(overrides.pilot, 'spacing')
    cfg.pilot.spacing = cfg.track_period / 8;
end
end

function merged = merge_struct(base, override)
merged = base;
fields = fieldnames(override);
for i = 1:numel(fields)
    name = fields{i};
    if isfield(base, name) && isstruct(base.(name)) && isstruct(override.(name))
        merged.(name) = merge_struct(base.(name), override.(name));
    else
        merged.(name) = override.(name);
    end
end
end
