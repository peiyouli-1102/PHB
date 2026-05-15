clear; clc; format long e;

% ===================== 积分 a =====================
% ∫_0^1 sin(x)/x dx,  x=0 处定义极限为 1
f_a = @(x) sinx_over_x(x);

% ===================== 积分 b =====================
% ∫_{-1}^1 (cos(x) - e^x)/sin(x) dx
% 处理 x=0 的奇点 (极限 -1), 避免减法精度损失
f_b = @(x) stable_b(x);

% ===================== 积分 c =====================
% 原题: ∫_0^∞ (x e^x)^{-1} dx
% 该积分在 x=0 发散, 按常见勘误此处应为 ∫_1^∞ (x e^x)^{-1} dx
% 作变量替换 x = 1/t → ∫_0^1 e^{-1/t}/t dt
% 若坚持原题, 可将区间改为 [1,∞) 并取消下面变换的注释
f_c = @(t) exp_m1_over_t(t);   % 在 [0,1] 上

% ========== 调用龙贝格算法, 计算 7 行 ==========
n_rows = 7;

fprintf('============= 积分 a: ∫_0^1 sin(x)/x dx =============\n');
R_a = romberg(f_a, 0, 1, n_rows);
print_romberg(R_a);

fprintf('\n============= 积分 b: ∫_{-1}^1 (cos x - e^x)/sin x dx =============\n');
R_b = romberg(f_b, -1, 1, n_rows);
print_romberg(R_b);

fprintf('\n============= 积分 c: ∫_1^∞ (x e^x)^{-1} dx  (变换后 ∫_0^1 e^{-1/t}/t dt) =============\n');
R_c = romberg(f_c, 0, 1, n_rows);
print_romberg(R_c);

% ==================== 辅助函数 ====================

function y = sinx_over_x(x)
    % 向量化处理, x=0 时返回 1
    y = ones(size(x));
    mask = (x ~= 0);
    y(mask) = sin(x(mask)) ./ x(mask);
end

function y = stable_b(x)
    % 稳定计算 (cos(x)-exp(x))/sin(x), 处理 x=0 奇点
    y = zeros(size(x));
    mask0 = (x == 0);
    y(mask0) = -1;                     % 极限值
    x_non = x(~mask0);
    % 使用恒等式: cos x - e^x = -2 sin^2(x/2) - (e^x - 1)
    %             sin x = 2 sin(x/2) cos(x/2)
    s = sin(x_non/2);
    c = cos(x_non/2);
    y(~mask0) = (-2*s.^2 - expm1(x_non)) ./ (2*s.*c);
end

function y = exp_m1_over_t(t)
    % 计算 e^{-1/t}/t, t=0 时极限为 0
    y = zeros(size(t));
    mask0 = (t == 0);
    y(mask0) = 0;
    t_non = t(~mask0);
    y(~mask0) = exp(-1./t_non) ./ t_non;
end

function print_romberg(R)
    % 打印下三角龙贝格阵列
    n = size(R,1);
    fprintf('龙贝格阵列 (行数 %d):\n', n);
    for i = 1:n
        for j = 1:i
            fprintf('%7.7f ', R(i,j));
            % fprintf('%7.10f ', R(i,j));
        end
        fprintf('\n');
    end
    fprintf('积分近似值 (对角线):  %.7f\n\n', R(n,n));
end

function R = romberg(f, a, b, n)
    % romberg - 龙贝格外推积分
    %   f   : 被积函数句柄 (支持向量输入)
    %   a,b : 积分区间
    %   n   : 要计算的阵列行数 (n×n 下三角)
    %   R   : n×n 龙贝格阵列 (只填充下三角, 上三角为0)
    
    R = zeros(n, n);
    h = b - a;
    % 第一行：1个区间的复化梯形值
    R(1,1) = (h/2) * (f(a) + f(b));
    
    % 逐行计算
    for i = 2:n
        h = h / 2;                     % 区间减半
        % 新增加的中点 (2^(i-2) 个)
        x_new = a + (2*(1:2^(i-2)) - 1) * h;
        sum_f = sum(f(x_new));
        % 更新梯形值
        R(i,1) = 0.5 * R(i-1,1) + h * sum_f;
        
        % 理查德森外推
        for j = 2:i
            R(i,j) = (4^(j-1) * R(i,j-1) - R(i-1,j-1)) / (4^(j-1) - 1);
        end
    end
end