clear; clc; format long;

% ---------- 1. 求解真值 (t=1) ----------
% 方程 x^2 - 1 + 2*exp(x) - 2*exp(-1) = 0 有两个根，正确根是 -1
f_true = @(x) x^2 - 1 + 2*exp(x) - 2*exp(-1);
options = optimset('Display','off','TolX',eps);
x_exact = fzero(f_true, -1, options);    % 初值改为 -1
fprintf('参考真解 x(1) = %.15f\n', x_exact);

% ---------- 2. 微分方程右端 ----------
f = @(t,x) (t - exp(-t)) ./ (x + exp(x));

% ---------- 3. 不同 N 循环 ----------
N_values = 2.^(3:8);          % N = 8,16,32,64,128,256
Errors = zeros(size(N_values));

for idx = 1:length(N_values)
    N = N_values(idx);
    h = 1 / N;
    t = linspace(0, 1, N+1)';
    x = zeros(N+1, 1);
    x(1) = 0;                 % 初值

    % --- 用经典 RK4 生成前 4 个初值 (x1, x2, x3) ---
    for i = 1:min(3, N)
        tn = t(i);
        xn = x(i);
        k1 = h * f(tn, xn);
        k2 = h * f(tn + h/2, xn + k1/2);
        k3 = h * f(tn + h/2, xn + k2/2);
        k4 = h * f(tn + h,   xn + k3);
        x(i+1) = xn + (k1 + 2*k2 + 2*k3 + k4)/6;
    end

    % --- Adams-Bashforth 4阶多步法 ---
    for i = 4:N
        fn   = f(t(i),   x(i));
        fn_1 = f(t(i-1), x(i-1));
        fn_2 = f(t(i-2), x(i-2));
        fn_3 = f(t(i-3), x(i-3));
        x(i+1) = x(i) + h/24 * (55*fn - 59*fn_1 + 37*fn_2 - 9*fn_3);
    end

    Errors(idx) = abs(x(end) - x_exact);
end

% ---------- 4. 计算收敛阶 ----------
Orders = NaN(size(Errors));
for idx = 2:length(N_values)
    Orders(idx) = log(Errors(idx-1) / Errors(idx)) / log(N_values(idx)/N_values(idx-1));
end

% ---------- 5. 输出表格 ----------
fprintf('\n   N          Error                Order\n');
fprintf('-----------------------------------------\n');
for idx = 1:length(N_values)
    if idx == 1
        fprintf('%4d    %.6e              -\n', N_values(idx), Errors(idx));
    else
        fprintf('%4d    %.6e        %.4f\n', N_values(idx), Errors(idx), Orders(idx));
    end
end