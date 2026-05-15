clear; close all;   % 清除工作区、命令窗口和图形
% clear; clc; close all;   % 清除工作区、命令窗口和图形

% 定义函数（匿名函数）
f = @(x) exp(x);

% 调用误差计算函数，同时获取线性样条和 clamped 三次样条的误差
[err_lin5, err_clamp5] = errcalculating(5, f);
[err_lin10, err_clamp10] = errcalculating(10, f);
[err_lin20, err_clamp20] = errcalculating(20, f);
[err_lin40, err_clamp40] = errcalculating(40, f);

% 计算收敛阶（基于线性样条）
ord_lin10 = log(err_lin5 / err_lin10) / log(2);
ord_lin20 = log(err_lin10 / err_lin20) / log(2);
ord_lin40 = log(err_lin20 / err_lin40) / log(2);

% 计算收敛阶（基于 clamped 样条）
ord_clamp10 = log(err_clamp5 / err_clamp10) / log(2);
ord_clamp20 = log(err_clamp10 / err_clamp20) / log(2);
ord_clamp40 = log(err_clamp20 / err_clamp40) / log(2);

% 输出表格
T = table([5;10;20;40], ...
          [err_lin5; err_lin10; err_lin20; err_lin40], ...
          [NaN ; ord_lin10; ord_lin20; ord_lin40], ...
          [err_clamp5; err_clamp10; err_clamp20; err_clamp40], ...
          [NaN ;ord_clamp10; ord_clamp20; ord_clamp40], ...
          'VariableNames', {'n', 'Linear_Error', 'Linear_Order', 'Clamped_Error', 'Clamped_Order'});
disp(T);

function [err1,err2] = errcalculating(n,f)
    for i = 1:n+1
        x_nodes(i) = (i-1)/(n);
    end
    y_nodes = f(x_nodes);  
    
    p1 = @(x) 0;
    p2 = @(x) 0;
    p3 = @(x) 0;
    
    % 误差计算
    xi = linspace((x_nodes(1)+x_nodes(2))/2,(x_nodes(n)+x_nodes(n+1))/2,n-1);
    yi_true = f(xi);        % 真实值
    
    yi_linear = hw3_linear_spline(x_nodes, y_nodes, xi);
    yi_clamped = hw3_clamped_cubic_spline(x_nodes, y_nodes, 1, exp(1), xi);
    
    err_linear = abs(yi_linear - yi_true);
    err_clamped = abs(yi_clamped - yi_true);
    
    err1 = max(err_linear);
    err2 = max(err_clamped);
end
