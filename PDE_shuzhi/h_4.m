% by cr0wchen
% for solving problem about textbook p143 7.2
clc; clear; close all;

f = @(x, y, t) sin(pi .* x) .* cos(pi .* y) .* exp(-(pi.^2) ./ 8 .* t);
ft0 = @(x, y)sin(pi .* x) .* cos(pi .* y);

J = 40;
N = 1600;
h = 1 / J;
tau = 1 / N; r = tau / (h^2);
x = 0 + [1:J - 1] * h; y = 0 + [1:J - 1] * h;
[X, Y] = meshgrid(x, y);

F = zeros(J - 1, J - 1, N + 1);
F(:, :, 1) = ft0(X, Y);
e = ones(J - 1, 1);
A = spdiags([e, -2 * e, e], [-1, 0, 1], J - 1, J - 1);
A = 1/16 * kron(A, spdiags([ones(J - 1, 1)], [0], J - 1, J - 1));
B = spdiags([e, -2 * e, e], [-1, 0, 1], J - 1, J - 1);
B(1, 1) = -1; B(J - 1, J - 1) = -1;
B = 1/16 * kron(spdiags([ones(J - 1, 1)], [0], J - 1, J - 1), B);
I = spdiags([ones(J - 1, 1)], [0], J - 1, J - 1);
I = kron(I, I);
temp = zeros((J - 1) * (J - 1), 1);

for n = 1:N
    temp = reshape(F(:, :, n), [], 1);
    temp = (I - r * A - r * B) \ temp;
    F(:, :, n + 1) = reshape(temp, [], J - 1);
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
