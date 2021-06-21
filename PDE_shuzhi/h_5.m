% for solving problem about p157
clc; clear; close all;

ut = @(x, t) cos(4 * pi .* t) .* sin(4 * pi .* x) + (sin(8 * pi .* t) .* sin(8 * pi .* x)) / (8 * pi);
ft0 = @(x) sin(4 * pi .* x);

J = 400; N = 500;
h = 1 / J; tao = 1 / N;
r = tao / h;
x = [1:J - 1] * h; t = [0:N] * tao;
U = zeros(J - 1, N + 1);
U(:, 1) = ft0(x);

e = ones(J - 1, 1);
U(:, 2) = tao * sin(8 * pi .* x') + sin(4 * pi .* x');

A = spdiags([e, -2 * e, e], [-1, 0, 1], J - 1, J - 1);

for n = 2:N
    U(:, n + 1) = r^2 * A * U(:, n) + 2 * U(:, n) - U(:, n - 1);
end

[T, X] = meshgrid(t, x);
figure
mesh(X, T, ut(X, T))
title("精确解")
figure
mesh(X, T, U);
title("数值解")
figure
mesh(X, T, U - ut(X, T));
title("误差")
