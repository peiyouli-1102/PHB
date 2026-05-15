function yi = hw3_linear_spline(x, y, xi)
% 线性样条插值
% 输入:
%   x - 节点横坐标 (递增)
%   y - 节点纵坐标
%   xi - 待插值点 (标量或向量)
% 输出:
%   yi - 插值结果
    n = length(x);
    n_xi = length(xi);
    yi = zeros(size(xi));
    for k = 1:n_xi
        xk = xi(k);
        if xk <= x(1)
            yi(k) = y(1);
        elseif xk >= x(end)
            yi(k) = y(end);
        else
            % 找到区间索引
            j = find(x > xk, 1, 'first');
            i = j - 1;
            t = (xk - x(i)) / (x(j) - x(i));
            yi(k) = (1 - t) * y(i) + t * y(j);
        end
    end
end