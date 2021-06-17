% for solving problem about p143 ADI
clc; clear; close all;

f = @(x, y, t) sin(pi .* x) .* cos(pi .* y) .* exp(-(pi.^2) ./ 8 .* t);
ft0 = @(x, y)sin(pi .* x) .* cos(pi .* y);

J = 40;
N = 1600;
h = 1 / J;
tao = 1 / N;
r = tao / (2 * h^2) / 16;
x = 0 + [1:J - 1] * h;
y = 0 + [1:J - 1] * h;
[X, Y] = meshgrid(x, y);

sz = J - 1
e = ones(sz, 1);
q_A = spdiags([-e * r, (1 + 2 * r) * e, -e * r], [-1, 0, 1], sz, sz);
q_B = spdiags([e * r, (1 - 2 * r) * e, e * r], [-1, 0, 1], sz, sz);
q_B(1, 1) = 1 - r;
q_B(sz, sz) = 1 - r;
h_A = q_A;
h_A(1, 1) = 1 + r;
h_A(sz, sz) = 1 + r;
h_B = q_B;

F = zeros(sz, sz, N + 1);
F(:, :, 1) = ft0(X, Y);
hf = zeros(sz, sz);
for n = 1:N
    hf = q_A' \ (F(:, :, n)' * q_B');
    F(:, :, n + 1) = h_A \ (hf' * h_B);
    
end

figure
mesh(X, Y, f(X, Y, 1))
title("精确解")
figure
mesh(X, Y, F(:, :, N + 1));
title("数值解")
figure
mesh(X, Y, F(:, :, N + 1) - f(X, Y, 1));
title("误差")
