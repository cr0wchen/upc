%by cr0wchen
% for solving problem about textbook p175 3 yin
clc; clear; close all;

J = 500; N = 1000;
h = 1 / J; tao = 1 / N;
a = -2;
r = a * tao / h;
x = [1:J - 1] * h; t = [0:N - 1] * tao;

U = zeros(J - 1, N);
U(:, 1) = sin(2 * pi * x) + 1;

e = ones(J - 1, 1);
A = spdiags([(1 - r) * e, r * e], [0, 1], J - 1, J - 1);

bouVal = zeros(J - 1, 1); bouVal(J - 1) = -r;
for n = 1:N - 1
    U(:, n + 1) = A \ (U(:, n) + bouVal);
end

[T, X] = meshgrid(t, x);

figure
mesh(X, T, U);
