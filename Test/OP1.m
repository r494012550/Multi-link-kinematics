% 多连杆机构向量法分析与轨迹优化（最终修正版）
%clear; clc; close all;

%% 1. 参数初始化
n_links = 4;                % 连杆数量
link_lengths = [1, 2, 1.5, 3];  % 初始连杆长度
initial_angles = [0, pi/4, pi/3];  % 初始关节角度（n_links-1个）
omega = 1;                  % 驱动角速度(rad/s)
sim_time = linspace(0, 2*pi, 100); % 仿真时间序列

% 优化参数边界约束
lb = [0.5*ones(1,n_links), -2*pi*ones(1,n_links-1)]; % 长度下限0.5m，角度±360°
ub = [5*ones(1,n_links), 2*pi*ones(1,n_links-1)];    % 长度上限5m

%% 2. 主优化流程
% 将长度和角度参数合并为优化变量
initial_params = [link_lengths, initial_angles];

% 配置优化选项
opt_options = optimoptions('fmincon',...
    'Display', 'iter',...     % 显示迭代过程
    'Algorithm', 'sqp',...    % 序列二次规划算法
    'MaxFunctionEvaluations', 1e4);

% 执行优化（添加非线性约束）
[optimal_params, fval] = fmincon(@trajectory_error, initial_params,...
    [], [], [], [], lb, ub, @closure_constraint, opt_options);

%% 3. 结果可视化
% 获取优化前后的轨迹数据
[init_x, init_y] = simulate_mechanism(link_lengths, initial_angles, sim_time, omega);
[opt_x, opt_y] = simulate_mechanism(optimal_params(1:n_links), optimal_params(n_links+1:end), sim_time, omega);

% 目标轨迹（示例：直线y=0.5x）
target_y = 0.5 * opt_x;

% 绘制轨迹对比
figure('Position', [100, 100, 800, 600])
plot(opt_x, opt_y, 'b-', 'LineWidth', 2, 'DisplayName','优化后轨迹');
hold on;
plot(opt_x, target_y, 'r--', 'LineWidth', 1.5, 'DisplayName','目标轨迹');
plot(init_x, init_y, 'g:', 'LineWidth', 1.5, 'DisplayName','初始轨迹');
title('多连杆机构轨迹优化对比');
xlabel('X位置 (m)'); ylabel('Y位置 (m)');
legend('Location','best'); grid on;
axis equal;

%% ------------------ 子函数定义 ------------------
% 目标函数：轨迹跟踪误差
function cost = trajectory_error(params)
    % 参数分解
    [lengths, angles] = parse_params(params);
    
    % 仿真获取末端轨迹
    [x, y] = simulate_mechanism(lengths, angles, evalin('base','sim_time'), evalin('base','omega'));
    
    % 计算目标轨迹误差（示例：直线y=0.5x）
    target_y = 0.5 * x;
    cost = sum((y - target_y).^2); % 最小二乘误差
end

% 非线性约束：机构闭合条件
function [c, ceq] = closure_constraint(params)
    [lengths, angles] = parse_params(params);
    
    % 计算各关节坐标
    [joint_x, joint_y] = calculate_joints(lengths, angles);
    
    % 闭合条件约束：首尾坐标重合
    ceq = [joint_x(end) - joint_x(1);  % X方向闭合
           joint_y(end) - joint_y(1)]; % Y方向闭合
    c = []; % 无不等式约束
end

% 机构运动学仿真（修正索引问题）
function [x_end, y_end] = simulate_mechanism(lengths, angles, t, omega)
    x_end = zeros(size(t));
    y_end = zeros(size(t));
    
    for i = 1:length(t)
        % 计算当前时刻各关节角度（假设第一个关节为驱动）
        current_angles = angles + omega * t(i);
        
        % 计算各节点坐标（修正循环逻辑）
        [x, y] = calculate_joints(lengths, current_angles);
        
        % 记录末端位置
        x_end(i) = x(end);
        y_end(i) = y(end);
    end
end

% 计算各关节坐标（核心运动学修正）
function [joint_x, joint_y] = calculate_joints(lengths, angles)
    n = length(lengths);
    joint_x = zeros(1, n+1); % 包含基座坐标
    joint_y = zeros(1, n+1);
    
    cumulative_angle = 0; % 初始角度
    
    % 修正循环逻辑：循环次数等于角度数量
    for k = 1:length(angles)
        cumulative_angle = cumulative_angle + angles(k);
        joint_x(k+1) = joint_x(k) + lengths(k)*cos(cumulative_angle);
        joint_y(k+1) = joint_y(k) + lengths(k)*sin(cumulative_angle);
    end
    
    % 处理最后一个连杆（无新角度）
    if n > length(angles)
        joint_x(end) = joint_x(end-1) + lengths(end)*cos(cumulative_angle);
        joint_y(end) = joint_y(end-1) + lengths(end)*sin(cumulative_angle);
    end
end

% 参数解析工具函数
function [lengths, angles] = parse_params(params)
    n_links = evalin('base', 'n_links');
    lengths = params(1:n_links);
    angles = params(n_links+1:end);
end