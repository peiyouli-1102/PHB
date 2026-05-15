function solve_ode_rkf45
% 使用 RKF45 自适应方法求解 y' = exp(y*x) + cos(y - x), y(1) = 3

% ---------- 初始条件与参数 ----------
x0 = 1; y0 = 3;
h = 0.01;               % 初值步长
tol = 1e-8;             % 局部误差容限
x_end_overflow = 1e10;  % 解溢出阈值
max_steps = 50000;      % 最大步数安全限制

% ---------- RKF45 系数 (Fehlberg 4(5)) ----------
c = [0; 1/4; 3/8; 12/13; 1; 1/2];
a = zeros(6,5);
a(2,1) = 1/4;
a(3,1) = 3/32;        a(3,2) = 9/32;
a(4,1) = 1932/2197;   a(4,2) = -7200/2197;  a(4,3) = 7296/2197;
a(5,1) = 439/216;     a(5,2) = -8;          a(5,3) = 3680/513;   
a(5,4) = -845/4104;
a(6,1) = -8/27;       a(6,2) = 2;           a(6,3) = -3544/2565;
a(6,4) = 1859/4104;   a(6,5) = -11/40;

% 4阶解系数 (低阶)
b4 = [25/216, 0, 1408/2565, 2197/4104, -1/5, 0];
% 5阶解系数 (高阶，用于推进)
b5 = [16/135, 0, 6656/12825, 28561/56430, -9/50, 2/55];

% ---------- 存储解 ----------
x_vals = x0;
y_vals = y0;

% ---------- 积分主循环 ----------
x = x0; y = y0;
for i = 1:max_steps
    if h < eps
        break;
    end
    
    % ---- 计算 6 个阶段 ----
    k = zeros(6,1);
    k(1) = f(x, y);
    for j = 2:6
        xj = x + c(j)*h;
        yj = y + h * (a(j,1:j-1) * k(1:j-1));
        k(j) = f(xj, yj);
    end
    
    % 4阶解与5阶解
    y4 = y + h * (b4 * k);
    y5 = y + h * (b5 * k);
    
    % 误差估计
    err = abs(y5 - y4);
    
    % ---------- 第二种步长策略 ----------
    if err <= tol
        % 接受步长
        x = x + h;
        y = y5;          % 使用高阶解推进
        x_vals = [x_vals, x];
        y_vals = [y_vals, y];
        
        % 溢出检测
        if abs(y) > x_end_overflow || isnan(y) || isinf(y)
            break;
        end
        
        % 计算新步长
        if err > 0
            h_new = 0.9 * h * (tol / err)^(1/5);
        else
            h_new = 2 * h;
        end
        % 限制步长变化范围
        h = min(5*h, max(0.2*h, h_new));
    else
        % 拒绝步长，根据误差缩小 h
        if err > 0
            h_new = 0.9 * h * (tol / err)^(1/5);
        else
            h_new = 0.5 * h;
        end
        h = max(0.2*h, h_new);  % 避免缩得太小
    end
    
    % 防止步长过小导致死循环
    if h < eps
        break;
    end
end

% ---------- 输出结果 ----------
fprintf('解的范围: [1, %.6f]\n', x_vals(end));

% 用户交互：线性插值
while true
    xq = input('请输入一个介于上述范围的值: ');
    if isempty(xq)
        break;
    end
    if xq < x_vals(1) || xq > x_vals(end)
        fprintf('输入超出范围，请重新输入.\n');
        continue;
    end
    % 两点线性插值
    yq = interp1(x_vals, y_vals, xq, 'linear');
    fprintf('插值结果 y(%.6f) = %.10f\n', xq, yq);
    break;
end

end

% ---------- 微分方程右端函数 ----------
function dydx = f(x, y)
    dydx = exp(y * x) + cos(y - x);
end