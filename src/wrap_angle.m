function theta = wrap_angle(theta)
%WRAP_ANGLE Wrap radians to [-pi, pi) without requiring a toolbox.
theta = mod(theta + pi, 2 * pi) - pi;
end
