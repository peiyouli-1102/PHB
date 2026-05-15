function yi = hw3_natural_cubic_spline(x, y, xi)
% 自然三次样条插值
% 输入输出同上
    n = length(x);
    h = diff(x);
    % 构建三对角矩阵（内部节点 M2...M_{n-1}）
    A = zeros(n-2, n-2);
    b = zeros(n-2, 1);
    for i = 2:n-1
        row = i-1;
        if i > 2
            A(row, row-1) = h(i-1);
        end
        A(row, row) = 2*(h(i-1) + h(i));
        if i < n-1
            A(row, row+1) = h(i);
        end
        b(row) = 6 * ( (y(i+1)-y(i))/h(i) - (y(i)-y(i-1))/h(i-1) );
    end
    M_int = A \ b;          % 内部节点的二阶导数
    M = zeros(n, 1);
    M(1) = 0;
    M(end) = 0;
    M(2:end-1) = M_int;

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