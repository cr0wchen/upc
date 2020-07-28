clc; clear; close all;

num = 40; %区间个数
n = 2; % 多项式的最高次数

%构造与多项式有关的A矩阵
mid = zeros(1, n + 1);
for i = 0:n
    mid(i + 1) = nchoosek(n, i);
end
A = diag(mid);
for i = 0:n
    for j = i + 1:n
        A(i + 1, j + 1) = (-1)^j * nchoosek(n, i)...
                         * nchoosek(n - i, n - i);
    end
end

%构造微分算子矩阵V
V = 1:n;
V = diag(V, -1);
% V = V(:, 1:n);

%构造与弱奇异积分有关的矩阵D
al = 1/2; %积分中的(t-s)的次数
D = zeros(1, n + 1);
s = 1;
for i = 1:n + 1
    s = s * (i - al);
    D(i) = factorial(i - 1) ./ s;
end
D = diag(D);

a = 0; b = 2; %左右两端点
t = linspace(a, b, num + 1); %等距节点
f = @(t) gamma(3) ./ gamma(2.5) * t.^1.5 -t.^2 - 1; %等式右端函数
la = 1 / gamma(0.5);
b = f(t)'; %方程组的常数项
MAT = zeros(length(t), n + 1); %c^T的系数矩阵
for i = 1:11
    MAT(i, :) = -A * power(t(i), 0:n)' + la * t(i)^(al) ...
                * A * V * D * power(t(i), 0:n)';
end
c = MAT \ b; %求解多项式系数
u = zeros(1, length(t));
for i = 1:length(t)
    u(i) = c' * A * power(t(i), 0:n)'; %得到数值解
end
ff = @(x) 1 + x.^2; %解析解的函数
plot(t, u, 'o')
hold on;
plot(t, ff(t), 'g')
title('微分方程数值解')
xlabel('t')
ylabel('u(t)')
legend('数值解', '真解');
wucha = abs(ff(t) - u);
wucha_per = sum(wucha) / length(wucha)%平均误差
