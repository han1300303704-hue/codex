function results = run_suite(cfg)
%RUN_SUITE Run all ablations, robustness scans, and write figures/data.

if nargin < 1
    cfg = default_config();
end
if ~exist(cfg.output_dir, 'dir')
    mkdir(cfg.output_dir);
end

scenarios = make_scenarios();
core = [];
for s = 1:numel(scenarios)
    fprintf('[%d/%d] Running %s with %d Monte Carlo trials...\n', ...
        s, numel(scenarios), scenarios(s).label, cfg.n_mc);
    result = run_monte_carlo(cfg, scenarios(s));
    if s == 1
        core = result;
    else
        core(s) = result;
    end
end

bit_scan = struct('bits', cfg.scan.phase_bits, 'loss_probability', zeros(size(cfg.scan.phase_bits)), ...
    'mean_gain', zeros(size(cfg.scan.phase_bits)));
for i = 1:numel(cfg.scan.phase_bits)
    scan_cfg = cfg;
    scan_cfg.n_mc = cfg.scan.n_mc;
    scan_cfg.hardware.phase_shifter_bits = cfg.scan.phase_bits(i);
    fprintf('[bit scan %d/%d] B=%d...\n', i, numel(cfg.scan.phase_bits), scan_cfg.hardware.phase_shifter_bits);
    scan_result = run_monte_carlo(scan_cfg, scenarios(end));
    bit_scan.loss_probability(i) = scan_result.loss_probability;
    bit_scan.mean_gain(i) = mean(scan_result.mean_gain);
end

pn_scan = struct('linewidth_hz', cfg.scan.phase_linewidth_hz, ...
    'loss_probability', zeros(size(cfg.scan.phase_linewidth_hz)), ...
    'mean_gain', zeros(size(cfg.scan.phase_linewidth_hz)));
for i = 1:numel(cfg.scan.phase_linewidth_hz)
    scan_cfg = cfg;
    scan_cfg.n_mc = cfg.scan.n_mc;
    scan_cfg.hardware.phase_noise_linewidth_hz = cfg.scan.phase_linewidth_hz(i);
    fprintf('[phase-noise scan %d/%d] linewidth=%g Hz...\n', i, numel(cfg.scan.phase_linewidth_hz), ...
        scan_cfg.hardware.phase_noise_linewidth_hz);
    scan_result = run_monte_carlo(scan_cfg, scenarios(end));
    pn_scan.loss_probability(i) = scan_result.loss_probability;
    pn_scan.mean_gain(i) = mean(scan_result.mean_gain);
end

results = struct('core', core, 'bit_scan', bit_scan, 'phase_noise_scan', pn_scan, ...
    'created_at', datestr(now, 30));
save(fullfile(cfg.output_dir, 'tracking_results.mat'), 'results', 'cfg');
save_result_figures(results, cfg);
end
