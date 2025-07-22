% 多连杆机构轨迹优化与运动动画
clear; clc; close all;

%% 1. 参数初始化
n_links = 4;                        % 连杆数量
link_lengths = [1, 2, 1.5, 3];      % 初始连杆长度(m)
initial_angles = [0, pi/4, pi/3];   % 初始关节角度(rad)
base_position = [0, 0];             % 基座位置(可修改)
omega = 1;                          % 驱动角速度(rad/s)
sim_time = linspace(0, 2*pi, 50);   % 仿真时间序列

%% 2. 机构运动学模型
% 获取初始位置数据
[init_x, init_y, init_joints] = simulate_mechanism(link_lengths, initial_angles, sim_time, omega, base_position);

%% 3. 创建动画窗口
figure('Position', [100 100 800 600])
hold on; axis equal; grid on;
xlabel('X位置 (m)'); ylabel('Y位置 (m)');
title('多连杆机构运动动画');

% 绘制目标轨迹
target_y = 0.5 * init_x;
plot(init_x, target_y, 'r--', 'LineWidth', 1.5, 'DisplayName','目标轨迹');

% 初始化动画元素
h_links = gobjects(1, n_links);     % 连杆图形对象
h_joints = gobjects(1, n_links+1);  % 关节标记
h_end = plot(NaN, NaN, 'b-', 'LineWidth', 2, 'DisplayName','末端轨迹'); % 末端轨迹

% 创建初始图形对象
for k = 1:n_links
    h_links(k) = line([0,0], [0,0], 'Color', '#0072BD', 'LineWidth', 2);
end
for k = 1:n_links+1
    h_joints(k) = plot(0, 0, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
end
legend('Location','northeast');

%% 4. 动画生成
for t_idx = 1:length(sim_time)
    % 获取当前时刻关节坐标
    current_joints = init_joints(:,:,t_idx);
    
    % 更新连杆位置
    for k = 1:n_links
        set(h_links(k),...
            'XData', [current_joints(1,k), current_joints(1,k+1)],...
            'YData', [current_joints(2,k), current_joints(2,k+1)]);
    end
    
    % 更新关节标记
    for k = 1:n_links+1
        set(h_joints(k),...
            'XData', current_joints(1,k),...
            'YData', current_joints(2,k));
    end
    
    % 更新末端轨迹
    new_x = [h_end.XData, current_joints(1,end)];
    new_y = [h_end.YData, current_joints(2,end)];
    set(h_end, 'XData', new_x, 'YData', new_y);
    
    % 设置动态视图范围
    axis([min(init_joints(1,:))-1, max(init_joints(1,:))+1,...
          min(init_joints(2,:))-1, max(init_joints(2,:))+1]);
    
    drawnow
    pause(0.05)  % 控制动画速度
end

%% ------------------ 核心函数 ------------------
function [x_end, y_end, all_joints] = simulate_mechanism(lengths, angles, t, omega, base)
    % 输入验证
    assert(length(lengths) == length(angles)+1,...
        '角度数量必须比连杆数量少1')
    
    % 预分配内存
    x_end = zeros(size(t));
    y_end = zeros(size(t));
    all_joints = zeros(2, length(lengths)+1, length(t)); % 3D数组存储所有关节位置
    
    for i = 1:length(t)
        % 计算当前角度
        current_angles = angles + omega * t(i);
        
        % 计算关节坐标
        [joints, ~] = calculate_joints(lengths, current_angles, base);
        
        % 记录数据
        x_end(i) = joints(1,end);
        y_end(i) = joints(2,end);
        all_joints(:,:,i) = joints;
    end
end

function [joints, cumulative_angles] = calculate_joints(lengths, angles, base)
    % 初始化存储数组
    n = length(lengths);
    joints = zeros(2, n+1);          % 二维坐标
    cumulative_angles = zeros(1, n); % 累计角度
    
    % 设置基座位置
    joints(:,1) = base';
    
    % 计算各关节位置
    current_angle = 0;
    for k = 1:n
        if k <= length(angles)
            current_angle = current_angle + angles(k);
        end
        cumulative_angles(k) = current_angle;
        
        joints(:,k+1) = joints(:,k) + lengths(k)*[cos(current_angle); 
                                                   sin(current_angle)];
    end
end