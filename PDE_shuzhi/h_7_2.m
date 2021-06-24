% by cr0wchen
% for solving problem about textbook p175 4(2)
clc; clear; close all;

J = 600; N = 100;
h = 3 / J; tau = 1 / N;
a = 1;
r = a * tau / h;
x = [1:J - 1] * h; t = [0:N - 1] * tau;

U = zeros(J - 1, N);
U(:, 1) = abs(x - 1);

e = ones(J - 1, 1);
A = spdiags([-r * e, (1 + r) * e], [-1, 0], J - 1, J - 1);

bouVal = zeros(J - 1, 1);
bouVal(1) = r;
for n = 1:N - 1
    U(:, n + 1) = A \ (U(:, n) + bouVal);
end

[T, X] = meshgrid(t, x);

figure
mesh(X, T, U);
