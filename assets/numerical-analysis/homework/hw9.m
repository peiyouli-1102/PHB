clear; clc; close all;

format short e

A = [-8003, 1999; 23988, -6004] ;
y0 = [1; 4] ;
h = 0.004 ;
N = 12 ;

% 显式欧拉方法（Explicit Euler / Forward Euler）

I_plus_A = eye(2) + h * A;
M = zeros(2, N+1);
M(:, 1) = y0;
y = [0; 0] ;
for i = 1:N
    y = I_plus_A * M(:, i) ;
    M(:, i+1) = y;
end

fprintf('显式欧拉方法：\n')
disp(M)

% 隐式欧拉方法（Implicit Euler / Backward Euler）

I_minus_A = eye(2) - h * A;
M = zeros(2, N+1);
M(:, 1) = y0;
for i = 1:N
    y = I_minus_A \ M(:, i);
    M(:, i+1) = y;
end

fprintf('隐式欧拉方法：\n')
disp(M)

% 梯形公式方法（Trapezoidal Rule，Crank–Nicolson ）

I_minus_A = eye(2) - (h/2) * A;
I_plus_A = eye(2) + (h/2) * A;
M = zeros(2, N+1);
M(:, 1) = y0;
for i = 1:N
    y = I_minus_A \ (I_plus_A * M(:, i));
    M(:, i+1) = y;
end

fprintf('梯形公式方法：\n')
disp(M)




