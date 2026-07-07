function channel = spherical_channel(r, theta, cfg)
%SPHERICAL_CHANNEL Near-field ULA line-of-sight channel and diagnostics.
% The ULA lies on the x-axis and broadside points along the positive y-axis.
% The default distance model uses the second-order Taylor/Fresnel
% approximation rho_n ~= r - x_n sin(theta) + x_n^2 cos^2(theta)/(2r).

r = max(real(r), cfg.channel.min_range);
n = (0:cfg.n_ant-1).' - (cfg.n_ant - 1) / 2;
x = n * cfg.spacing;
user_x = r * sin(theta);
user_y = r * cos(theta);
rho_exact = sqrt((x - user_x).^2 + user_y.^2);
rho = element_distances_taylor(r, theta, x, cfg, rho_exact);
k = 2 * pi / cfg.lambda;

% Remove a common propagation phase so the array response is well scaled.
steering = exp(-1j * k * (rho - r)) ./ sqrt(cfg.n_ant);
amplitude = cfg.channel.path_gain_at_1m ./ rho;
h = amplitude .* exp(-1j * k * (rho - r));

channel = struct('positions', x, 'distances', rho, 'exact_distances', rho_exact, ...
    'distance_error', rho - rho_exact, 'distance_model', cfg.channel.distance_model, ...
    'steering', steering, 'h', h, 'user_position', [user_x; user_y]);
end

function rho = element_distances_taylor(r, theta, x, cfg, rho_exact)
model = 'exact';
if isfield(cfg, 'channel') && isfield(cfg.channel, 'distance_model')
    model = lower(cfg.channel.distance_model);
end
if strcmp(model, 'exact')
    rho = rho_exact;
    return;
end

s = sin(theta);
c = cos(theta);
rho = r - x * s + (x .^ 2) * (c ^ 2) / (2 * r);
rho = max(real(rho), eps);
end
