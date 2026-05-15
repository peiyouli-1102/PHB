clear; close all;   % 清除工作区、命令窗口和图形
% clear; clc; close all;   % 清除工作区、命令窗口和图形

% 定义函数（匿名函数）
f = @(x) 1./(1+x.^2);

n = input('N=');

for i = 1:n+1
    % 计算第 i 个点的 x 坐标
    xnodes1(i) = 5 - 10 * (i-1)/(n);
    ynodes1(i) = f(xnodes1(i));
    xnodes2(i) = -5 * cos((2*i-1)/(2*n+2)*pi);
    ynodes2(i) = f(xnodes2(i));
end

p1 = @(x) 0;
p2 = @(x) 0;

for i = 1:n+1
    % 计算点 x_i 对应基函数
    p_t1 = @(x) 1;
    p_t2 = @(x) 1;
    for j = 1:n+1
        if j ~= i
            p_t1 = @(x) p_t1(x) .* ( x - xnodes1(j)) / ( xnodes1(i) - xnodes1(j) );
            p_t2 = @(x) p_t2(x) .* ( x - xnodes2(j)) / ( xnodes2(i) - xnodes2(j) );
        end
    end
    p1 = @(x) p1(x) + ynodes1(i) * p_t1(x) ;
    p2 = @(x) p2(x) + ynodes2(i) * p_t2(x) ;
end

serr1 = 0;
serr2 = 0;
for i = 1:101
    x_t = (i-1) / 10 - 5; 
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

x_plot = linspace(-5, 5, 500);
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