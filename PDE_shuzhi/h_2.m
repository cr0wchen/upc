clc; clear; close all;

a = -1;
b = 1;
c = -2;
d = 2;
N = 1000;

h1 = (b - a) / N;
h2 = (d - c) / N;
uf = @(x, y) sin(pi .* (x + 1) ./ 2) .* sin(pi .* (y + 2) ./ 2);%真解
f = @(x, y) (pi.^2 ./ 2) .* sin(pi .* (x + 1) ./ 2) .* sin(pi .* (y + 2) ./ 2);

x_h = (a + h1 * [1:N - 1]);
y_h = (c + h2 * [1:N - 1]);

% [x, y] = meshgrid(x_h, y_h);
% x = reshape(x, [], 1);
% y = reshape(y, [], 1);

x = kron(x_h', ones(N - 1, 1));
y = kron(ones(N - 1, 1), y_h');

S = (N - 1) * (N - 1);
eL = ones(S, 1);
A1 = spdiags([eL, -2 * eL, eL], [-(N - 1), 0, N - 1], S, S);
A1 = -1 * A1 / (h1^2);
A2 = spdiags([ones(N - 1, 1), -2 * ones(N - 1, 1), ones(N - 1, 1)], [-1, 0, 1], N - 1, N - 1);
A2 = kron(spdiags([ones(N - 1, 1)], [0], N - 1, N - 1), A2);

A2 = -1 * A2 / (h2^2);
u = (A1 + A2) \ f(x, y);
[x_m, y_m] = meshgrid(x_h, y_h);
u_t = uf(x_m, y_m);
u_Z = reshape(u, [], N - 1);

figure;
mesh(x_m, y_m, uf(x_m, y_m));
title("真解")
figure;
mesh(x_m, y_m, reshape(u_t, [], N - 1));
title("数值解");
figure;
mesh(x_m, y_m, u_Z - u_t)
title("误差")
