clear; close all;   % 清除工作区、命令窗口和图形
% clear; clc; close all;   % 清除工作区、命令窗口和图形

% 定义被积函数与积分区间
f1 = @(x) sin(x);      a1 = 0; b1 = 4;      exact1 = 1 - cos(4);   % 真值
f2 = @(x) sin(x);      a2 = 0; b2 = 2*pi;   exact2 = 0;            % 真值

% 准备 N = 2^k, k = 1,...,12
kmax = 12;
Nvals = 2.^(1:kmax);   % 行向量

% 预存结果
Trap_I1 = zeros(1, kmax);   Trap_I2 = zeros(1, kmax);
Simp_I1 = zeros(1, kmax);   Simp_I2 = zeros(1, kmax);
Trap_err1 = zeros(1, kmax); Trap_err2 = zeros(1, kmax);
Simp_err1 = zeros(1, kmax); Simp_err2 = zeros(1, kmax);

% 循环计算
for idx = 1:kmax
    N = Nvals(idx);
    % 梯形公式
    Trap_I1(idx) = CompTrapezoidal(f1, a1, b1, N);
    Trap_I2(idx) = CompTrapezoidal(f2, a2, b2, N);
    % Simpson 公式
    Simp_I1(idx) = CompSimpson(f1, a1, b1, N);
    Simp_I2(idx) = CompSimpson(f2, a2, b2, N);
    % 误差 (取绝对值)
    Trap_err1(idx) = abs(Trap_I1(idx) - exact1);
    Trap_err2(idx) = abs(Trap_I2(idx) - exact2);
    Simp_err1(idx) = abs(Simp_I1(idx) - exact1);
    Simp_err2(idx) = abs(Simp_I2(idx) - exact2);
end

% 计算收敛阶 (从 k=2 开始)
Ord_Trap1 = zeros(1, kmax-1);  Ord_Simp1 = zeros(1, kmax-1);
Ord_Trap2 = zeros(1, kmax-1);  Ord_Simp2 = zeros(1, kmax-1);

for k = 2:kmax
    N_old = Nvals(k-1);  N_now = Nvals(k);
    err_old = Trap_err1(k-1);  err_now = Trap_err1(k);
    Ord_Trap1(k-1) = log(err_old/err_now) / log(N_now/N_old);
    
    err_old = Simp_err1(k-1);  err_now = Simp_err1(k);
    Ord_Simp1(k-1) = log(err_old/err_now) / log(N_now/N_old);
    
    err_old = Trap_err2(k-1);  err_now = Trap_err2(k);
    Ord_Trap2(k-1) = log(err_old/err_now) / log(N_now/N_old);
    
    err_old = Simp_err2(k-1);  err_now = Simp_err2(k);
    Ord_Simp2(k-1) = log(err_old/err_now) / log(N_now/N_old);
end

% 显示结果表 
fprintf('\n========== 积分 I1 = ∫_0^4 sin(x) dx，真值 = %.12f ==========\n', exact1);
fprintf('  k      N       梯形近似值       梯形误差      Simpson近似值      Simpson误差\n');
for k = 1:kmax
    fprintf('%3d  %6d   %.12f   %.4e   %.12f   %.4e\n', ...
        k, Nvals(k), Trap_I1(k), Trap_err1(k), Simp_I1(k), Simp_err1(k));
end

fprintf('\n收敛阶 (I1):\n');
fprintf('  k    N_old -> N_new   梯形阶    Simpson阶\n');
for k = 2:kmax
    fprintf('%3d  %6d -> %6d   %6.3f    %6.3f\n', ...
        k, Nvals(k-1), Nvals(k), Ord_Trap1(k-1), Ord_Simp1(k-1));
end

fprintf('\n========== 积分 I2 = ∫_0^2pi sin(x) dx，真值 = %.13g ==========\n', exact2);
fprintf('  k      N       梯形近似值       梯形误差      Simpson近似值      Simpson误差\n');
for k = 1:kmax
    fprintf('%3d  %6d   %.13g   %.4e   %.13g   %.4e\n', ...
        k, Nvals(k), Trap_I2(k), Trap_err2(k), Simp_I2(k), Simp_err2(k));
end

fprintf('\n收敛阶 (I2):\n');
fprintf('  k    N_old -> N_new   梯形阶    Simpson阶\n');
for k = 2:kmax
    fprintf('%3d  %6d -> %6d   %6.3f    %6.3f\n', ...
        k, Nvals(k-1), Nvals(k), Ord_Trap2(k-1), Ord_Simp2(k-1));
end

% 可绘制双对数误差图
figure;
loglog(Nvals, Trap_err1, 'o-', 'LineWidth', 1.5, 'DisplayName', '梯形法 I1');
hold on;
loglog(Nvals, Simp_err1, 's-', 'LineWidth', 1.5, 'DisplayName', 'Simpson法 I1');
loglog(Nvals, Trap_err2, '^-', 'LineWidth', 1.5, 'DisplayName', '梯形法 I2');
loglog(Nvals, Simp_err2, 'd-', 'LineWidth', 1.5, 'DisplayName', 'Simpson法 I2');
xlabel('N (子区间数)'); ylabel('绝对误差');
legend('Location', 'southwest'); grid on;
title('复化求积公式的误差收敛行为');

function I = CompTrapezoidal(f, a, b, N)
% 复化梯形公式
% f - 被积函数
% a, b - 积分区间
% N - 子区间个数
    x = linspace(a, b, N+1);
    y = f(x);
    h = (b - a) / N;
    
    % 显式计算内部点的和 (y(2) 到 y(N))
    sum_inner = 0;
    for i = 2:N
        sum_inner = sum_inner + y(i);
    end
    
    I = h * (sum_inner + (y(1) + y(N+1))/2);
end

function I = CompSimpson(f, a, b, N)
% 复化 Simpson 公式
% f - 被积函数
% a, b - 积分区间
% N - 子区间个数 (必须为偶数)
    if mod(N, 2) ~= 0
        error('N 必须为偶数');
    end
    x = linspace(a, b, N+1);
    y = f(x);
    h = (b - a) / N;
    
    % 计算奇数下标点（不包括首尾）: y(3), y(5), ..., y(N-1)
    sum_odd = 0;
    for i = 3:2:N-1
        sum_odd = sum_odd + y(i);
    end
    
    % 计算偶数下标点（不包括首尾）: y(2), y(4), ..., y(N)
    sum_even = 0;
    for i = 2:2:N
        sum_even = sum_even + y(i);
    end
    
    I = h/3 * (y(1) + 4*sum_even + 2*sum_odd + y(N+1));
end

