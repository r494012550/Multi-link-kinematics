clc; clear;

% 已知常数
X = 2130; Y = 350;
L1 = 280; L2 = 1588.5; L3 = 940;
L4 = 840; L5 = 850; L6 = 845;

% θ1 范围与初始化
theta1_deg = 0:1:360;
theta1_rad = deg2rad(theta1_deg);
n = length(theta1_rad);

% 输出初始化
S_list = nan(1, n);
sol_store = nan(n, 6);  % [theta2, theta3, theta4, theta5, theta6, S]

% 初始猜测
x0 = [0.5, -0.5, 0.2, 0.2, pi, 1000];  % [θ2,θ3,θ4,θ5,θ6,s]

% 设置 fsolve 选项
opts = optimoptions('fsolve', 'Display', 'off', 'FunctionTolerance', 1e-8);

% 逐点求解
for k = 1:n
    theta1 = theta1_rad(k);
    
    fun = @(x) real_imag_eqns(theta1, x, L1, L2, L3, L4, L5, L6, X, Y);
    
    if k > 1
        x0 = sol_store(k-1, :);  % 用前一解做初始猜测
    end

    [x_sol, ~, exitflag] = fsolve(fun, x0, opts);

    if exitflag > 0
        sol_store(k, :) = x_sol;
        S_list(k) = x_sol(6);
    else
        fprintf("解失败: θ1 = %.1f°\n", theta1_deg(k));
        x0 = [0.5, -0.5, 0.2, 0.2, pi, 1000];  % 重置猜测防止错误传播
    end
end

% 绘图
plot(theta1_deg, S_list, 'b-', 'LineWidth', 1.5);
xlabel('\theta_1 (°)');
ylabel('滑块位移 S (mm)');
title('S vs θ_1');
grid on;
