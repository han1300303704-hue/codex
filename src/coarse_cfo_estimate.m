function cfo = coarse_cfo_estimate(y, pilot_spacing)
%COARSE_CFO_ESTIMATE Common CFO estimate from adjacent spatial pilot vectors.

if size(y, 2) < 2
    cfo = 0;
    return;
end
correlation = sum(conj(y(:, 1:end-1)) .* y(:, 2:end), 1);
cfo = angle(sum(correlation)) / (2 * pi * pilot_spacing);
end
