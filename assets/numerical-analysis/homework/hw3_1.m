clear; close all;   % 清除工作区、命令窗口和图形
% clear; clc; close all;   % 清除工作区、命令窗口和图形

% 定义函数（匿名函数）
f = @(x) 1./(1+ x.^2);

x_nodes = [0, 5/3, 10/3, 5];
y_nodes = f(x_nodes);  

p1 = @(x) 0;
p2 = @(x) 0;
p3 = @(x) 0;

% 精细插值点（用于绘图和误差计算）
xi = linspace(min(x_nodes), max(x_nodes), 500);
yi_true = f(xi);        % 真实值

% 调用三种插值函数
yi_linear = hw3_linear_spline(x_nodes, y_nodes, xi);
yi_natural = hw3_natural_cubic_spline(x_nodes, y_nodes, xi);
yi_hermite = hw3_hermite_spline(x_nodes, y_nodes, xi);

% 计算误差（绝对值）
err_linear = abs(yi_linear - yi_true);
err_natural = abs(yi_natural - yi_true);
err_hermite = abs(yi_hermite - yi_true);

fprintf("Max Error of linear spline       : %e\n", max(err_linear));
fprintf("Max Error of natural cubic spline: %e\n", max(err_natural));
fprintf("Max Error of hermite spline      : %e\n", max(err_hermite));

x_plot = linspace(0, 5, 500);
y_f = f(x_plot);
y_p1 = p1(x_plot);
y_p2 = p2(x_plot);

% 绘图：插值曲线对比
figure('Position', [100, 100, 1200, 500]);

subplot(1,2,1);
plot(xi, yi_true, 'k-', 'LineWidth', 2, 'DisplayName', '真实函数 f(x)');
hold on;
plot(x_nodes, y_nodes, 'ro', 'MarkerSize', 8, 'LineWidth', 1.5, 'DisplayName', '节点');
plot(xi, yi_linear, 'b-', 'LineWidth', 1.2, 'DisplayName', '线性样条');
plot(xi, yi_natural, 'g-', 'LineWidth', 1.2, 'DisplayName', '自然三次样条');
plot(xi, yi_hermite, 'm-', 'LineWidth', 1.2, 'DisplayName', '三次 Hermite');
legend('Location', 'best');
xlabel('x'); ylabel('y');
title('插值曲线对比');
grid on;

subplot(1,2,2);
semilogy(xi, err_linear, 'b-', 'LineWidth', 1.2, 'DisplayName', '线性样条误差');
hold on;
semilogy(xi, err_natural, 'g-', 'LineWidth', 1.2, 'DisplayName', '自然三次样条误差');
semilogy(xi, err_hermite, 'm-', 'LineWidth', 1.2, 'DisplayName', '三次 Hermite 误差');
legend('Location', 'best');
xlabel('x'); ylabel('绝对误差');
title('误差曲线（对数坐标）');
grid on;
