clear; close all;   % 清除工作区、命令窗口和图形
% clear; clc; close all;   % 清除工作区、命令窗口和图形

% 定义函数（匿名函数）
f = @(x) 1./(1+ 25.*x.^2);

n = input('N=');

for i = 1:n+1
    % 计算第 i 个点的 x 坐标
    xnodes1(i) = 1 - 2 * (i-1)/(n);
    ynodes1(i) = f(xnodes1(i));
    xnodes2(i) = -cos((2*i-1)/(2*n+2)*pi);
    ynodes2(i) = f(xnodes2(i));
end

p1 = @(x) 0;
p2 = @(x) 0;

% 初始化 n×n 的零矩阵 存储差商
D1 = zeros(n+1);
D2 = zeros(n+1);

% 计算差商
for i = 1:n+1
    D1(i,1) = ynodes1(i); 
    D2(i,1) = ynodes2(i);
    p_t1 = @(x) 1;
    p_t2 = @(x) 1;
    for j = 2:i
        D1(i,j) = (D1(i,j-1)-D1(i-1,j-1))/(xnodes1(i)-xnodes1(i-j+1)) ;
        D2(i,j) = (D2(i,j-1)-D2(i-1,j-1))/(xnodes2(i)-xnodes2(i-j+1)) ;
        p_t1 = @(x) p_t1(x) .* ( x - xnodes1(j-1)) ;
        p_t2 = @(x) p_t2(x) .* ( x - xnodes2(j-1)) ;
    end
    p1 = @(x) p1(x) + D1(i,i) .* p_t1(x) ;
    p2 = @(x) p2(x) + D2(i,i) .* p_t2(x) ;
end

serr1 = 0;
serr2 = 0;
for i = 1:101
    x_t = (i-1) / 50 - 1; 
    err1(i) = abs(f(x_t)-p1(x_t));
    err2(i) = abs(f(x_t)-p2(x_t));
    if err1(i) > serr1
        serr1 = err1(i);
    end
    if err2(i) > serr2
        serr2 = err2(i);
    end
end

fprintf("Max Error of grid (1) : %e\n", serr1);
fprintf("Max Error of grid (2) : %e\n", serr2);

x_plot = linspace(-1, 1, 500);
y_f = f(x_plot);
y_p1 = p1(x_plot);
y_p2 = p2(x_plot);

figure;
plot(x_plot, y_f, 'b-', 'LineWidth', 1.5); hold on;
plot(x_plot, y_p1, 'g-', 'LineWidth', 1.5);hold on;
plot(x_plot, y_p2, 'r-', 'LineWidth', 1.5);
plot(xnodes1, ynodes1, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
plot(xnodes2, ynodes2, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
xlabel('x'); ylabel('y');
legend('原函数 f(x)', '插值多项式 P1(x)','插值多项式 P2(x)', '插值节点');
title(sprintf('%d 次拉格朗日插值', n));
grid on;