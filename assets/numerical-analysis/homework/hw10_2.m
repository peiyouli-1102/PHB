clear; clc; format long;

% ---------- 初值（扰动） ----------
x0 = 1e-4;                % 不再使用 0

% ---------- 微分方程右端 ----------
f = @(t,x) (t - exp(-t)) ./ (x + exp(x));

% ---------- 高精度参考解 ----------
options_ode = odeset('RelTol',1e-13, 'AbsTol',1e-13);
[~, x_ref] = ode45(f, [0 1], x0, options_ode);
x_exact = x_ref(end);
fprintf('高精度参考解 x(1) = %.15f\n\n', x_exact);

% ---------- 不同 N 的误差计算 ----------
N_values = 2.^(3:8);          % 8,16,32,64,128,256
Errors = zeros(size(N_values));

for idx = 1:length(N_values)
    N = N_values(idx);
    h = 1 / N;
    t = linspace(0, 1, N+1)';
    x = zeros(N+1, 1);
    x(1) = x0;               % 扰动后的初值

    % --- RK4 起步（计算 x1, x2, x3） ---
    for i = 1:min(3, N)
        tn = t(i);
        xn = x(i);
        k1 = h * f(tn, xn);
        k2 = h * f(tn + h/2, xn + k1/2);
        k3 = h * f(tn + h/2, xn + k2/2);
        k4 = h * f(tn + h,   xn + k3);
        x(i+1) = xn + (k1 + 2*k2 + 2*k3 + k4)/6;
    end

    % --- Adams-Bashforth 4 阶多步法 ---
    for i = 4:N
        fn   = f(t(i),   x(i));
        fn_1 = f(t(i-1), x(i-1));
        fn_2 = f(t(i-2), x(i-2));
        fn_3 = f(t(i-3), x(i-3));
        x(i+1) = x(i) + h/24 * (55*fn - 59*fn_1 + 37*fn_2 - 9*fn_3);
    end

    Errors(idx) = abs(x(end) - x_exact);
end

% ---------- 收敛阶 ----------
Orders = NaN(size(Errors));
for idx = 2:length(N_values)
    Orders(idx) = log(Errors(idx-1) / Errors(idx)) / ...
                  log(N_values(idx) / N_values(idx-1));
end

% ---------- 输出表格 ----------
fprintf('   N          Error                Order\n');
fprintf('-----------------------------------------\n');
for idx = 1:length(N_values)
    if idx == 1
        fprintf('%4d    %.6e              -\n', N_values(idx), Errors(idx));
    else
        fprintf('%4d    %.6e        %.4f\n', ...
                N_values(idx), Errors(idx), Orders(idx));
    end
end