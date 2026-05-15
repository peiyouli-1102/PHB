clear; clc; format long;

% ---------- 被积函数定义 ----------
f1 = @(x) exp(-x.^2);
f2 = @(x) 1./(1 + x.^2);
f3 = @(x) 1./(2 + cos(x));

% ---------- 积分区间与精确值 ----------
a1=0; b1=1;   I1_exact = sqrt(pi)/2 * erf(1);       % 0.746824132812427
a2=0; b2=4;   I2_exact = atan(4);                   % 1.32581766366803
a3=0; b3=2*pi; I3_exact = 2*pi/sqrt(3);             % 3.62759872846844

% 被积函数与区间存为 cell，方便循环
funcs = {f1, f2, f3};
exact = [I1_exact, I2_exact, I3_exact];
intervals = [a1, b1; a2, b2; a3, b3];

% ---------- 参数设置 ----------
k_vals = 1:7;
N_vals = 2.^k_vals;            % N = 2,4,8,16,32,64,128
nN = length(N_vals);

% 预存误差
err_trap = zeros(nN, 3);
err_gauss = zeros(nN, 3);

% ---------- 主循环 ----------
for idx = 1:nN
    N = N_vals(idx);
    for s = 1:3          % s=1,2,3 分别对应 I1, I2, I3
        f = funcs{s};
        a = intervals(s,1);
        b = intervals(s,2);

        % 复化梯形
        trap_val = composite_trapezoidal(f, a, b, N);
        err_trap(idx, s) = abs(trap_val - exact(s));

        % 复化 3 点 Gauss
        gauss_val = composite_gauss3(f, a, b, N);
        err_gauss(idx, s) = abs(gauss_val - exact(s));
    end
end

% ---------- 计算收敛阶 ----------
% 阶的定义：ln(Error_old / Error_now) / ln(N_now / N_old)
% 因为 N 每次翻倍，分母 ln(2)，实际阶 = log2(Error_old / Error_now)
order_trap = zeros(nN, 3);
order_gauss = zeros(nN, 3);
for s = 1:3
    order_trap(2:end, s) = log2(err_trap(1:end-1, s) ./ err_trap(2:end, s));
    order_gauss(2:end, s) = log2(err_gauss(1:end-1, s) ./ err_gauss(2:end, s));
end

% ---------- 表格输出 ----------
fprintf('==================== 复化梯形公式 ====================\n');
print_table(N_vals, err_trap, order_trap);

fprintf('\n================== 复化 3 点 Gauss 公式 ==================\n');
print_table(N_vals, err_gauss, order_gauss);

% ========== 辅助函数 ==========

function T = composite_trapezoidal(f, a, b, N)
    h = (b - a) / N;
    x = linspace(a, b, N+1);
    y = f(x);
    T = h * (0.5*y(1) + sum(y(2:end-1)) + 0.5*y(end));
end

function G = composite_gauss3(f, a, b, N)
    % 3 点 Gauss-Legendre 节点与权重 (区间 [-1,1])
    xg = [-sqrt(3/5), 0, sqrt(3/5)];
    wg = [5/9, 8/9, 5/9];
    h = (b - a) / N;
    mids = a + h/2 : h : b - h/2;   % 各子区间中点
    % 构造所有积分节点 (向量化)
    nodes = mids + xg' * (h/2);     % 3×N 矩阵
    fvals = f(nodes);               % 函数值 3×N
    G = sum(wg * fvals) * (h/2);    % wg:1×3, fvals:3×N → 1×N 求和
end

function print_table(N_vals, err, order)
    fprintf('%-8s %-18s %-8s %-18s %-8s %-18s %-8s\n', ...
        'N', 'I1误差', '阶', 'I2误差', '阶', 'I3误差', '阶');
    fprintf('-------------------------------------------------------------\n');
    for i = 1:length(N_vals)
        if i == 1
            fprintf('%-8d %-18.4e %-8s %-18.4e %-8s %-18.4e %-8s\n', ...
                N_vals(i), err(i,1), '-', err(i,2), '-', err(i,3), '-');
        else
            fprintf('%-8d %-18.4e %-8.3f %-18.4e %-8.3f %-18.4e %-8.3f\n', ...
                N_vals(i), err(i,1), order(i,1), ...
                err(i,2), order(i,2), err(i,3), order(i,3));
        end
    end
end