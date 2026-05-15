function yi = hw3_hermite_spline(x, y, xi)
% 三次 Hermite 插值（导数由中心差分估算）
    n = length(x);
    % 估算导数
    dy = zeros(1, n);
    for i = 1:n
        if i == 1
            dy(i) = (y(2) - y(1)) / (x(2) - x(1));
        elseif i == n
            dy(i) = (y(n) - y(n-1)) / (x(n) - x(n-1));
        else
            dy(i) = (y(i+1) - y(i-1)) / (x(i+1) - x(i-1));
        end
    end

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
            hi = x(j) - x(i);
            t = (xk - x(i)) / hi;
            H0 = 2*t^3 - 3*t^2 + 1;
            H1 = -2*t^3 + 3*t^2;
            H2 = t^3 - 2*t^2 + t;
            H3 = t^3 - t^2;
            yi(k) = y(i)*H0 + y(j)*H1 + dy(i)*hi*H2 + dy(j)*hi*H3;
        end
    end
end