%RUN_EXPERIMENT Main entry point for the near-field predictive beam tracker.
% Results are written to results/ as .mat, .fig and .png files.

root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'src'));
cfg = default_config();
cfg.output_dir = fullfile(root, 'results');
results = run_suite(cfg); %#ok<NASGU>
fprintf('Completed. Results are available in %s\n', cfg.output_dir);
