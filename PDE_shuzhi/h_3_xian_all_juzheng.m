% for solving problem about p141
clc; clear; close all;
ut = @(t, x) exp(-pi^2 .* t) .* cos(pi * x) + (1 - cos(t));
f = @(t) sin(t);

h = 1/40;
tao = 1/3200;
r = 1/2;

t = 0:tao:1;
x = 0:h:1;
J = length(x);
N = length(t);

F = zeros(J, N);
F(:, 1) = cos(pi * x);
e = ones(J, 1);
A = spdiags([r * e, (1 - 2 * r) * e, r * e], [-1, 0, 1], J, J);

for n = 1:N -1
    F(:, n + 1) = A * F(:, n) + tao * f(t(n));
    F(1, n + 1) = ut(t(n + 1), 0);
    F(J, n + 1) = ut(t(n + 1), 1);
end

[t_m, x_m] = meshgrid(t, x);
u_t = ut(t_m, x_m);
figure;
mesh(t_m, x_m, u_t)
title("真解")
figure;
mesh(t_m, x_m, F)
title("数值解")
figure;
mesh(t_m, x_m, abs(u_t - F))
% mesh(t_m, x_m, u_t - F)
title("误差")
