clc; clear; close all;

J = 100; N = 500;
h = 1 / J; tao = 1 / N;
a = -2;
r = a * tao / h;
x = [1:J - 1] * h; t = [1:N - 1] * tao;
U = zeros(J + 1, N + 1);
U(:, 1) = sin(2 * [0, x, 1]') + 1;

e = ones(J + 1, 1);
A = spdiags([e, -e], [0, 1], J + 1, J + 1);

for n = 1:N
    % U(J + 1, n) = 1;
    U(:, n + 1) = r * A * U(:, n) + U(:, n);
end

size(U)
[T, X] = meshgrid([0, t, 1], [0, x, 1]);

figure
mesh(X, T, U);
