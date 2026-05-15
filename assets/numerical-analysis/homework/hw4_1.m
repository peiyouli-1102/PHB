clear; close all;   % 清除工作区、命令窗口和图形

% 定义函数句柄
u = @(x) 1 ./ (5 - 4*cos(x));
v = @(x) sin(x/2);

% 密集网格（用于计算精确误差，画图用）
x_fine = linspace(0, 2*pi, 2000)';
x_fine = x_fine(2:end-1);   % 去掉 0,2π 点，避免重复

N_list = [4, 8, 16, 32, 64];
colors = lines(length(N_list));

% 选择要处理的函数（可改为 u 或 v）
f = v;
fname = 'v';

% 存储误差（最大绝对误差）
err_proj = zeros(size(N_list));
err_interp = zeros(size(N_list));

% 绘图准备
figure('Position', [100 100 1200 500]);

% 子图1：投影误差
subplot(1,2,1);
hold on;
for i = 1:length(N_list)
    N = N_list(i);
    P_f = FourierProjection(f, N, x_fine);
    err = abs(f(x_fine) - P_f);
    err_proj(i) = max(err);          % 存储最大误差
    plot(x_fine, err, 'Color', colors(i,:), 'LineWidth', 1.5);
end
hold off;
xlim([0 2*pi]);
set(gca, 'YScale', 'log');
xlabel('x'); ylabel('|f(x) - P_N f(x)|');
title(sprintf('Projection error for %s', fname));
legend(arrayfun(@(n) sprintf('N=%d', n), N_list, 'UniformOutput', false));
grid on;

% 子图2：插值误差
subplot(1,2,2);
hold on;
for i = 1:length(N_list)
    N = N_list(i);
    I_f = FourierInterpolation(f, N, x_fine);
    err = abs(f(x_fine) - I_f);
    err_interp(i) = max(err);        % 存储最大误差
    plot(x_fine, err, 'Color', colors(i,:), 'LineWidth', 1.5);
end
hold off;
xlim([0 2*pi]);
set(gca, 'YScale', 'log');
xlabel('x'); ylabel('|f(x) - I_N f(x)|');
title(sprintf('Interpolation error for %s', fname));
legend(arrayfun(@(n) sprintf('N=%d', n), N_list, 'UniformOutput', false));
grid on;

% 第二个图形：投影和插值函数本身
figure('Position', [100 100 1200 500]);

subplot(1,2,1);
hold on;
for i = 1:length(N_list)
    N = N_list(i);
    P_f = FourierProjection(f, N, x_fine);
    plot(x_fine, P_f, 'Color', colors(i,:), 'LineWidth', 1.5);
end
hold off;
xlim([0 2*pi]);
xlabel('x'); ylabel('P_N f(x)');
title(sprintf('Projection for %s', fname));
legend(arrayfun(@(n) sprintf('N=%d', n), N_list, 'UniformOutput', false));
grid on;

subplot(1,2,2);
hold on;
for i = 1:length(N_list)
    N = N_list(i);
    I_f = FourierInterpolation(f, N, x_fine);
    plot(x_fine, I_f, 'Color', colors(i,:), 'LineWidth', 1.5);
end
hold off;
xlim([0 2*pi]);
xlabel('x'); ylabel('I_N f(x)');
title(sprintf('Interpolation for %s', fname));
legend(arrayfun(@(n) sprintf('N=%d', n), N_list, 'UniformOutput', false));
grid on;

% ---------- 计算收敛阶并输出表格 ----------
% 投影误差的收敛阶
order_proj = NaN(1, length(N_list));
for i = 2:length(N_list)
    order_proj(i) = log(err_proj(i-1) / err_proj(i)) / log(2);
end

% 插值误差的收敛阶
order_interp = NaN(1, length(N_list));
for i = 2:length(N_list)
    order_interp(i) = log(err_interp(i-1) / err_interp(i)) / log(2);
end

% 输出表格
T = table(N_list', err_proj', order_proj', err_interp', order_interp', ...
    'VariableNames', {'N', 'Proj_Error', 'Proj_Order', 'Interp_Error', 'Interp_Order'});
disp(' ');
disp(['收敛阶分析 (函数: ', fname, ')']);
disp(T);

% ---------- 以下为原始函数定义，保持不变 ----------
function P_f = FourierProjection(f, N, x_fine)
    % 手搓版：不使用 fft，直接计算傅里叶系数
    K = N/2;                    % 最大频率
    M = 2^nextpow2(10*max(N,64));  % 采样点数（足够密）
    x_sample = linspace(0, 2*pi, M+1)';
    x_sample = x_sample(1:end-1);
    f_sample = f(x_sample);
    dx = 2*pi / M;
    
    % 计算傅里叶系数 hat{f}_k (k = -K .. K)
    coeffs = zeros(1, 2*K+1);
    for idx = 1:2*K+1
        k = idx - K - 1;        % k 从 -K 到 K
        sum_val = 0;
        for j = 1:M
            sum_val = sum_val + f_sample(j) * exp(-1i * k * x_sample(j));
        end
        coeffs(idx) = sum_val / M;
    end
    
    % 在 x_fine 上重构
    P_f = zeros(size(x_fine));
    for idx = 1:2*K+1
        k = idx - K - 1;
        P_f = P_f + coeffs(idx) * exp(1i * k * x_fine);
    end
    P_f = real(P_f);
end

function I_f = FourierInterpolation(f, N, x_fine)
    L = N + 1;
    x_nodes = (0:N) * 2*pi / L;
    f_nodes = f(x_nodes);
    K = N/2;
    
    % 计算插值系数 d_k (k = -K .. K)
    d = zeros(1, 2*K+1);
    for idx = 1:2*K+1
        k = idx - K - 1;
        sum_val = 0;
        for j = 1:L
            sum_val = sum_val + f_nodes(j) * exp(-1i * k * x_nodes(j));
        end
        d(idx) = sum_val / L;
    end
    
    % 重构
    I_f = zeros(size(x_fine));
    for idx = 1:2*K+1
        k = idx - K - 1;
        I_f = I_f + d(idx) * exp(1i * k * x_fine);
    end
    I_f = real(I_f);
end