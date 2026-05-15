function yi = hw3_clamped_cubic_spline(x, y, d1, dn, xi)
% 完备三次样条（固定边界导数）
% 输入：
%   x  - 节点横坐标（递增）
%   y  - 节点纵坐标
%   d1 - 左端点一阶导数 S'(x1)
%   dn - 右端点一阶导数 S'(xn)
%   xi - 待插值点（标量或向量）
% 输出：
%   yi - 插值结果
    n = length(x);
    h = diff(x);
    
    % 构建三对角矩阵 A (n×n) 和右端项 b (n×1)
    A = zeros(n, n);
    b = zeros(n, 1);
    
    % 内部节点方程 (i = 2,...,n-1)
    for i = 2:n-1
        A(i, i-1) = h(i-1);
        A(i, i)   = 2*(h(i-1) + h(i));
        A(i, i+1) = h(i);
        b(i) = 6 * ( (y(i+1)-y(i))/h(i) - (y(i)-y(i-1))/h(i-1) );
    end
    
    % 左端点边界条件
    A(1, 1) = 2*h(1);
    A(1, 2) = h(1);
    b(1) = 6 * ( (y(2)-y(1))/h(1) - d1 );
    
    % 右端点边界条件
    A(n, n-1) = h(n-1);
    A(n, n)   = 2*h(n-1);
    b(n) = 6 * ( dn - (y(n)-y(n-1))/h(n-1) );
    
    % 求解 M = S''(x_i)
    M = A \ b;
    
    % 计算插值
    n_xi = length(xi);
    yi = zeros(size(xi));
    for k = 1:n_xi
        xk = xi(k);
        if xk <= x(1)
            yi(k) = y(1);
        elseif xk >= x(end)
            yi(k) = y(end);
        else
            j = find(x > xk, 1, 'first');
            i = j - 1;
            hi = h(i);
            Mi = M(i);
            Mj = M(j);
            term1 = ( (x(j)-xk)^3 * Mi + (xk-x(i))^3 * Mj ) / (6*hi);
            term2 = ( y(i) - Mi*hi^2/6 ) * (x(j)-xk)/hi;
            term3 = ( y(j) - Mj*hi^2/6 ) * (xk-x(i))/hi;
            yi(k) = term1 + term2 + term3;
        end
    end
end