% 机构参数初始化（需要根据实际机构尺寸修改）
L1 = 237;      % 曲柄长度
L2 = 1150;      % 连杆长度
L3 = 1020;    % 延长杆长度
L4 = 800;      % 摇杆长度
L6 = 1000;      % 滑块导轨长度
x0 = 1250;      % 原点x坐标
y0 = 350;      % 原点y坐标

% 运动参数初始化
w1 = -pi/6;  % 曲柄角速度(rad/s)
t = 0:0.01:12; % 时间序列
dt = 0.01;    % 时间步长

% 预分配存储数组
s = zeros(size(t));

for i = 1:length(t)
    % ========== 位移分析 ==========
    theta1 = w1*t(i) + 3*pi/2; % 曲柄角度
    
    % 计算辅助参数
    A = L1*cos(theta1) - x0;
    B = L1*sin(theta1) + y0;
    C = (A^2 + B^2 + L4^2 - L2^2)/(2*L4);
    
    % 计算theta4（注意处理复数解）
    discriminant = A^2 + B^2 - C^2;
    if discriminant < 0
        error('机构无法装配，请检查参数。t=%.2f', t(i));
    end
    theta4 = 2*atan((B - sqrt(discriminant))/(A - C));
    
    % 计算theta2
    theta2 = atan2(B + L4*sin(theta4), A + L4*cos(theta4));
    
    % 计算其他角度
    theta3 = theta4 + 53*pi/45;
    numerator = L1*cos(theta1) + L3*cos(theta3);
    if abs(numerator/L6) > 1
        error('超出acos定义域，t=%.2f', t(i));
    end
    theta6 = acos(numerator/L6);
    
    % 计算滑块位移
    s(i) = L1 + L3 + L6 + L1*sin(theta1) + L3*sin(theta3) - L6*sin(theta6);
end

% ========== 速度加速度计算 ==========
% 数值微分法
v = gradient(s, dt);      % 一阶导数求速度
a = gradient(v, dt);      % 二阶导数求加速度

% 绘制结果
figure;
subplot(3,1,1)
plot(t, s)
title('滑块位移')
xlabel('时间(s)')
ylabel('位移(m)')

subplot(3,1,2)
plot(t, v)
title('滑块速度')
xlabel('时间(s)')
ylabel('速度(m/s)')

subplot(3,1,3)
plot(t, a)
title('滑块加速度')
xlabel('时间(s)')
ylabel('加速度(m/s²)')