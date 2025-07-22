% 机构参数
L1 = 208; L2 = 1588.5; L3 = 940; L4 = 840; L5 = 850; L6 = 845;
X0 = 2130; Y0 = 350;

% 配置求解器选项
options = optimoptions('fsolve', 'Display', 'off', 'FunctionTolerance', 1e-6);

% 定义角度范围和步长 (高分辨率)
theta1_deg = 8.7:1:360+8.7;        % 1度步长 (361个点)

theta1_rad = deg2rad(theta1_deg); % 角度转换为弧度
n = length(theta1_rad);

% 预分配内存
s_values = zeros(1, n);
theta_sol = zeros(n, 6);      % 存储所有变量解 [θ2,θ3,θ4,θ5,θ6,s]

% 初始猜测 (基于机构特性)
x0 = [0.5, -0.5, 0, 0, pi, 1000];  % [θ2,θ3,θ4,θ5,θ6,s]

% 循环求解每个θ1对应的s
for i = 1:n
    % 定义当前θ1的方程组
    eqns = @(x) [
        L1*cos(theta1_rad(i)) + L2*cos(x(1)) + L3*cos(x(2)) - X0;
        L1*sin(theta1_rad(i)) + L2*sin(x(1)) + L3*sin(x(2)) + Y0;
        L4*cos(x(3)) + L5*cos(x(4)) - L3*cos(x(2));
        L4*sin(x(3)) + L5*sin(x(4)) - L3*sin(x(2));
        L4*cos(x(3)) + L6*cos(x(5));
        L4*sin(x(3)) + L6*sin(x(5)) + x(6);
    ];
    
    % 使用前一点解作为初始猜测 (加速收敛)
    if i > 1
        x0 = theta_sol(i-1, :);
    end
    
    % 求解方程组 
    [x_sol, ~, exitflag] = fsolve(eqns, x0, options);
    
    % 存储有效解
    if exitflag > 0
        theta_sol(i, :) = x_sol;
        s_values(i) = x_sol(6);
    else
        error('求解失败于 θ1 = %.1f°', theta1_deg(i));
    end
end

% 创建行程函数 s(θ1) 使用三次样条插值
s_func = @(th1_deg) interp1(theta1_deg, s_values, th1_deg, 'spline');

% 计算速度 v(θ1) = ds/dθ1 (单位: mm/rad)
% 修正: 使用数值微分计算速度
dtheta_rad = diff(theta1_rad);
ds = diff(s_values);
v_values = ds ./ dtheta_rad;  % 速度值 (每个区间)

% 创建速度函数 (插值)
v_func = @(th1_deg) interp1(theta1_deg(1:end-1), v_values, th1_deg, 'spline', 'extrap');

% 计算加速度 a(θ1) = dv/dθ1 (单位: mm/rad²)
% 修正: 使用数值微分计算加速度
dv = diff(v_values);
a_values = dv ./ dtheta_rad(1:end-1);  % 加速度值

% 创建加速度函数 (插值)
a_func = @(th1_deg) interp1(theta1_deg(1:end-2), a_values, th1_deg, 'spline', 'extrap');

% ===== 结果验证和可视化 =====
% 计算关键位置的运动参数
test_angles = [0, 8.27 ,90 , 150 , 180, 188.27 , 210 , 270, 360];
fprintf('θ1(°)   s(mm)    v(mm/rad)   a(mm/rad²)\n');
fprintf('---------------------------------------\n');

for th = test_angles
    s_val = s_func(th);
    v_val = v_func(th);
    a_val = a_func(th);
    fprintf('%5.1f   %7.1f   %9.2f   %9.2f\n', th, s_val, v_val, a_val);
end

% 绘制运动曲线
figure('Position', [100, 100, 1200, 800]);

% 位移曲线
subplot(3,1,1);
plot(theta1_deg, s_values, 'b-', 'LineWidth', 1.5);
xlabel('θ_1 (°)');
ylabel('位移 s (mm)');
title('滑块行程曲线');
grid on;

% 速度曲线
subplot(3,1,2);
v_plot = v_func(theta1_deg(1:end-1));  % 有效范围
plot(theta1_deg(1:end-1), v_plot, 'r-', 'LineWidth', 1.5);
xlabel('θ_1 (°)');
ylabel('速度 v (mm/rad)');
title('滑块速度曲线');
grid on;

% 加速度曲线
subplot(3,1,3);
a_plot = a_func(theta1_deg(1:end-2));  % 有效范围
plot(theta1_deg(1:end-2), a_plot, 'g-', 'LineWidth', 1.5);
xlabel('θ_1 (°)');
ylabel('加速度 a (mm/rad²)');
title('滑块加速度曲线');
grid on;

% ===== 实用函数 =====
% 获取完整行程数据
function [s, v, a] = get_motion_data(theta1_deg)
    s = s_func(theta1_deg);
    v = v_func(theta1_deg);
    a = a_func(theta1_deg);
end

% 转换为时间导数 (输入转速 deg/s)
function [v_t, a_t] = convert_to_time_derivatives(theta1_deg, dtheta1_dt_deg)
    dtheta1_dt_rad = deg2rad(dtheta1_dt_deg);  % 转换为 rad/s
    
    % 获取位置导数 (rad)
    v_rad = v_func(theta1_deg);
    a_rad = a_func(theta1_deg);
    
    % 转换为时间导数
    v_t = v_rad * dtheta1_dt_rad;           % mm/s
    a_t = a_rad * (dtheta1_dt_rad)^2;       % mm/s²
end