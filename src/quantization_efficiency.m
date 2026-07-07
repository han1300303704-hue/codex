function eta = quantization_efficiency(bits)
%QUANTIZATION_EFFICIENCY Approximate phase-only quantization gain ceiling.
%   ETA is the large-array coherent power efficiency relative to continuous
%   phase control.  It is useful for reporting hardware-normalized gain:
%       G_hw = G_continuous_normalized / ETA.
%
%   For B-bit uniform phase quantization,
%       ETA = (sin(pi/2^B)/(pi/2^B))^2.

if isinf(bits) || isnan(bits)
    eta = 1;
    return;
end
levels = 2^bits;
x = pi / levels;
eta = (sin(x) / x)^2;
end
