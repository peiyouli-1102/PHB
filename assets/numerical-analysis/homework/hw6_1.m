clear; close all;   % 清除工作区、命令窗口和图形
% clear; clc; close all;   % 清除工作区、命令窗口和图形

% 定义函数（匿名函数）

u = @(x) log(x);
v = @(x) tan(x);
w = @(x) sin(x^2 + 1/3*x);

RichardsonExtrapolation(u, 3, 3, 1);
RichardsonExtrapolation(v, asin(0.8), 4, 1);
RichardsonExtrapolation(w, 0, 5, 1);

function f_1 = RichardsonExtrapolation(f, x, M, h)
    D = zeros(M+1);
    for k = 1:(M+1)
        step = h / (2 ^ k);
        D(k,1) = (f(x + step) - f(x - step)) / (2 * step);
    end
    pow4 = 4.^(1:M+1);   % 预计算4的幂次
    for m = 2:(M+1)
        for k = m:(M+1)
            D(k,m) = D(k,m-1) + (D(k,m-1) - D(k-1,m-1)) / (pow4(m-1) - 1);
        end
    end
    for i = 1:M+1
        for j = 1:i
            fprintf('%8.4f ', D(i,j));   % 每个数字固定宽度，自然对齐
        end
        fprintf('\n');   % 换行
    end
    f_1 = D(M+1,M+1);
end

