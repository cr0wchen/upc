clc; clear; close all;

num = 40; %�������
n = 2; % ����ʽ����ߴ���

%���������ʽ�йص�A����
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

%����΢�����Ӿ���V
V = 1:n;
V = diag(V, -1);
% V = V(:, 1:n);

%����������������йصľ���D
al = 1/2; %�����е�(t-s)�Ĵ���
D = zeros(1, n + 1);
s = 1;
for i = 1:n + 1
    s = s * (i - al);
    D(i) = factorial(i - 1) ./ s;
end
D = diag(D);

a = 0; b = 2; %�������˵�
t = linspace(a, b, num + 1); %�Ⱦ�ڵ�
f = @(t) gamma(3) ./ gamma(2.5) * t.^1.5 -t.^2 - 1; %��ʽ�Ҷ˺���
la = 1 / gamma(0.5);
b = f(t)'; %������ĳ�����
MAT = zeros(length(t), n + 1); %c^T��ϵ������
for i = 1:11
    MAT(i, :) = -A * power(t(i), 0:n)' + la * t(i)^(al) ...
                * A * V * D * power(t(i), 0:n)';
end
c = MAT \ b; %������ʽϵ��
u = zeros(1, length(t));
for i = 1:length(t)
    u(i) = c' * A * power(t(i), 0:n)'; %�õ���ֵ��
end
ff = @(x) 1 + x.^2; %������ĺ���
plot(t, u, 'o')
hold on;
plot(t, ff(t), 'g')
title('΢�ַ�����ֵ��')
xlabel('t')
ylabel('u(t)')
legend('��ֵ��', '���');
wucha = abs(ff(t) - u);
wucha_per = sum(wucha) / length(wucha)%ƽ�����
