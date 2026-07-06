function save_result_figures(results, cfg)
%SAVE_RESULT_FIGURES Save the required curves as MATLAB figures and PNGs.

labels = {results.core.label};
time_ms = results.core(1).time_s * 1e3;
colors = lines(numel(results.core));

fig = figure('Visible', cfg.figure_visible, 'Color', 'w', 'Name', 'Beam gain');
hold on;
for i = 1:numel(results.core)
    plot(time_ms, results.core(i).mean_gain, 'LineWidth', 1.4, 'Color', colors(i, :));
end
grid on; ylim([0 1.05]); xlabel('时间 (ms)'); ylabel('归一化波束增益');
legend(labels, 'Location', 'best'); title('近场预测波束跟踪增益');
save_figure_pair(fig, cfg.output_dir, 'beam_gain'); close(fig);

fig = figure('Visible', cfg.figure_visible, 'Color', 'w', 'Name', 'Tracking error');
tiledlayout(3, 1, 'TileSpacing', 'compact');
nexttile; hold on;
for i = 1:numel(results.core)
    plot(time_ms, results.core(i).position_rmse_m, 'LineWidth', 1.3, 'Color', colors(i, :));
end
grid on; ylabel('位置 RMSE (m)'); title('跟踪误差');
nexttile; hold on;
for i = 1:numel(results.core)
    plot(time_ms, results.core(i).range_rmse_m, 'LineWidth', 1.3, 'Color', colors(i, :));
end
grid on; ylabel('距离 RMSE (m)');
nexttile; hold on;
for i = 1:numel(results.core)
    plot(time_ms, rad2deg(results.core(i).angle_rmse_rad), 'LineWidth', 1.3, 'Color', colors(i, :));
end
grid on; xlabel('时间 (ms)'); ylabel('角度 RMSE (deg)');
legend(labels, 'Location', 'best');
save_figure_pair(fig, cfg.output_dir, 'tracking_error'); close(fig);

fig = figure('Visible', cfg.figure_visible, 'Color', 'w', 'Name', 'Lock loss');
losses = [results.core.loss_probability];
slot_outages = [results.core.slot_outage_probability];
bar([losses(:), slot_outages(:)]); grid on; ylim([0 1]);
set(gca, 'XTick', 1:numel(labels), 'XTickLabel', labels, 'XTickLabelRotation', 20);
ylabel('概率'); title('固定参数下的失锁概率');
legend({'曾经失锁概率', '时隙失锁占比'}, 'Location', 'best');
save_figure_pair(fig, cfg.output_dir, 'lock_loss_probability'); close(fig);

fig = figure('Visible', cfg.figure_visible, 'Color', 'w', 'Name', 'Robustness scans');
tiledlayout(1, 2, 'TileSpacing', 'compact');
nexttile;
plot(results.bit_scan.bits, results.bit_scan.loss_probability, '-o', 'LineWidth', 1.5);
grid on; ylim([0 1]); xlabel('移相器量化位宽 B'); ylabel('失锁概率'); title('量化位宽扫描');
nexttile;
semilogx(results.phase_noise_scan.linewidth_hz, results.phase_noise_scan.loss_probability, '-s', 'LineWidth', 1.5);
grid on; ylim([0 1]); xlabel('相位噪声线宽 (Hz)'); ylabel('失锁概率'); title('相位噪声扫描');
save_figure_pair(fig, cfg.output_dir, 'robustness_scans'); close(fig);
end

function save_figure_pair(fig, output_dir, stem)
savefig(fig, fullfile(output_dir, [stem '.fig']));
saveas(fig, fullfile(output_dir, [stem '.png']));
end
