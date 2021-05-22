% for solving two-point boundary problem -d(xdu/dx)dx+x^2(du/dx)-xu=(x^2-1)cosx on [0,pi]
% with u(0)=u(pi)=0
% true solution u=sinx
clc; clear; close all;

a = 0;
b = pi;
N = 1000;
h = (b - a) / N;
x = a + [1:N - 1] * h;
fx = @(x) (x.^2 - 1) .* cos(x);
dx = [x + h / 2]';
ds = [x - h / 2]';
A = spdiags([dx, -(dx + ds), ds], [-1, 0, 1], N - 1, N - 1);
A = -1 * A / (h^2);
e = ones(N - 1, 1);
B = spdiags([-e, zeros(N - 1, 1), e], [-1, 0, 1], N - 1, N - 1);
B = spdiags([x.^2]', 0, N - 1, N - 1) * B / (2 * h);
C = spdiags([-x]', 0, N - 1, N - 1);
vector_f = fx(x)';
u = (A + B + C) \ vector_f;
% plot(x, u)

% hold on;
% plot(x, sin(x), '--r')

plot(x,u'-sin(x))
