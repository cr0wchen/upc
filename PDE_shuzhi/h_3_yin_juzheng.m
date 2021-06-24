% by cr0wchen
% for solving problem about textbook p141 7.1
clc; clear; close all;
ut = @(t, x) exp(-pi^2 .* t) .* cos(pi * x) + (1 - cos(t));
f = @(t) sin(t);

J = 40;
N = 1600;
h = 1 / J;
tau = 1 / N;
a = 1;
r = a * tau / (h^2);

t = [0:N - 1] * tau;
x = [1:J - 1] * h;

U = zeros(J - 1, N);
U(:, 1) = cos(pi * x);
e = ones(J - 1, 1);
A = spdiags([-r * e, (1 + 2 * r) * e, -r * e], [-1, 0, 1], J - 1, J - 1);
tVec = zeros(J - 1, 1);

for n = 1:N -1
    tVec(1) = r * ut(t(n), 0); tVec(J - 1) = r * ut(t(n), 1);
    U(:, n + 1) = A \ (U(:, n) + tau * f(t(n)) + tVec);
end

[t_m, x_m] = meshgrid(t, [0, x, 1]);
u_t = ut(t_m, x_m);
figure;
mesh(t_m, x_m, u_t)
title("真解")
figure;
mesh(t_m, x_m, [ut(t, 0); U; ut(t, 1)])
title("数值解")
figure;
% mesh(t_m, x_m, abs(u_t - U))
mesh(t_m, x_m, u_t - [ut(t, 0); U; ut(t, 1)])
title("误差")
