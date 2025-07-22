% 曲柄滑块机构运动分析
% 参数定义
r = 160;        % 曲柄半径 (mm)
l = 1050;       % 连杆长度 (mm)
ab = 800;      %连杆质心到滑块的距离 (mm)
omega = 2*pi;   % 曲柄角速度 (rad/s)  实际角速度 可能为一个函数 omega(t)
dt = 0.001;     % 时间步长 (s)

% 时间向量 (0 到 1 秒，完整周期)
t = 0:dt:1;     % 1秒对应完整周期(ω=2π rad/s)

% 曲柄角度 (随时间变化)
theta = omega * t + pi;

%连杆与曲柄的交点
xa = @(th)r*sin(th);
ya = @(th)r*cos(th);

% 滑块位置函数
yb = @(th) r*cos(th) + sqrt(l^2 - (r*sin(th)).^2);   %核心公式 论文上面人工推导的公式，应采用向量方程，由程序求解

%连杆质心位移与角位移
xc = xa(theta) * ab / l;
yc = ya(theta) + ( yb(theta) - ya(theta)) * (l - ab) / l;
theta2 = asin( xa(theta) / l);

%求导
v = gradient(yb(theta), dt); %速度 ----结果已验证合理
a = gradient(v, dt);  %加速度 ----结果存疑

%连杆的线速度
vxc = gradient(xc, dt);
vyc = gradient(yc, dt);
vsc = sqrt(vxc.^2 + vyc.^2);

%连杆的角速度
vtheta2 = gradient(theta2, dt);

% 显示关键时间点的运动参数
fprintf('时间(s)\t角度(deg)\t位置(mm)\t速度(mm/s)\t加速度(mm/s²)\t连杆线速度(mm/s)\t连杆角速度(rad/s)\n');
fprintf('------------------------------------------------------------------------------------------------------\n');
for ti = [0, 0.25, 0.463, 0.5, 0.75, 1]
    index = round(ti/dt) + 1;
    fprintf('%.3f\t%.2f\t\t%.2f\t\t%.2f\t\t%.2f\t\t%.2f\t\t%.2f\n', ...
            t(index), theta(index)/pi*180, yb(theta(index)), ...
            v(index), a(index), vsc(index), vtheta2(index));
end


% 分别绘制位置、速度和加速度
figure;
subplot(5,1,1);
plot(t, yb(theta));
title('滑块位置随时间变化');
xlabel('时间 (s)');
ylabel('位置 (mm)');

subplot(5,1,2);
plot(t, v);
title('滑块速度随时间变化');
xlabel('时间 (s)');
ylabel('速度 (mm/s)');

subplot(5,1,3); 
plot(t, a);
title('滑块加速度随时间变化');
xlabel('时间 (s)');
ylabel('加速度 (mm/s²)');

subplot(5,1,4); 
plot(t, vsc);
title('连杆线速度随时间变化');
xlabel('时间 (s)');
ylabel('线速度 (mm/s)');

subplot(5,1,5); 
plot(t, vtheta2);
title('连杆角速度随时间变化');
xlabel('时间 (s)');
ylabel('角速度 (rad/s)');