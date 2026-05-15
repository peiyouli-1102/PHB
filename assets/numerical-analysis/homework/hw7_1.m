clear; close all;

N_vals = [5,10,15,20,25,30,35,40];
fprintf("Exact : %e \n", 2/5 * atan(5));
for n = N_vals
    interp_int_approx(n);
end

function err = interp_int_approx(n)
    f = @(x) 1./(1 + 25*x.^2);
    I_exact = 2/5 * atan(5);

    % 生成节点
    xnodes1 = zeros(1, n+1);
    ynodes1 = zeros(1, n+1);
    xnodes2 = zeros(1, n+1);
    ynodes2 = zeros(1, n+1);
    for i = 1:n+1
        xnodes1(i) = 1 - 2*(i-1)/n;
        ynodes1(i) = f(xnodes1(i));
        xnodes2(i) = -cos( i / (n+2) * pi );   
        ynodes2(i) = f(xnodes2(i));
    end

    % 初始化多项式系数
    p1_coeff = zeros(1, n+1);
    p2_coeff = zeros(1, n+1);

    for i = 1:n+1
        L1 = 1;
        L2 = 1;
        for j = 1:n+1
            if j ~= i
                L1 = conv(L1, [1, -xnodes1(j)]);
                L2 = conv(L2, [1, -xnodes2(j)]);
            end
        end
        den1 = prod(xnodes1(i) - xnodes1([1:i-1, i+1:end]));
        den2 = prod(xnodes2(i) - xnodes2([1:i-1, i+1:end]));
        L1 = L1 / den1;
        L2 = L2 / den2;
        p1_coeff = p1_coeff + ynodes1(i) * L1;
        p2_coeff = p2_coeff + ynodes2(i) * L2;
    end

    I1 = polyval(polyint(p1_coeff), 1) - polyval(polyint(p1_coeff), -1);
    I2 = polyval(polyint(p2_coeff), 1) - polyval(polyint(p2_coeff), -1);
    err1 = abs(I1 - I_exact);
    err2 = abs(I2 - I_exact);

    fprintf("N = %d ; I1 = %.3e ; I2 = %.3e ; err1 = %.3e ; err2 = %.3e\n", ...
             n, I1, I2, err1, err2);
    err = err1;   % 只返回第一组误差（可根据需要修改）
end

