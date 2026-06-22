function channel = spherical_channel(r, theta, cfg)
%SPHERICAL_CHANNEL Near-field ULA line-of-sight channel and diagnostics.
% The ULA lies on the x-axis and broadside points along the positive y-axis.

r = max(real(r), cfg.channel.min_range);
n = (0:cfg.n_ant-1).' - (cfg.n_ant - 1) / 2;
x = n * cfg.spacing;
user_x = r * sin(theta);
user_y = r * cos(theta);
rho = sqrt((x - user_x).^2 + user_y.^2);
k = 2 * pi / cfg.lambda;

% Remove a common propagation phase so the array response is well scaled.
steering = exp(-1j * k * (rho - r)) ./ sqrt(cfg.n_ant);
amplitude = cfg.channel.path_gain_at_1m ./ rho;
h = amplitude .* exp(-1j * k * (rho - r));

channel = struct('positions', x, 'distances', rho, 'steering', steering, ...
    'h', h, 'user_position', [user_x; user_y]);
end
