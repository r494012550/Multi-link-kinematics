% 参数设置
r = 50;       % 曲柄半径(mm)
L = 150;      % 连杆长度(mm)
e = -30;      % 向左偏心(mm)
omega = 2*pi; % 角速度(rad/s, 对应60RPM)
theta = -10:0.1:370;   % 扩展角度范围观察死点行为

% 转换为弧度
theta_rad = deg2rad(theta);

% ===== 速度计算 =====
% 无偏心(e=0)
dsdtheta0 = -r*sin(theta_rad) - (r^2*sin(theta_rad).*cos(theta_rad))...
            ./sqrt(L^2 - (r*sin(theta_rad)).^2 + eps);
v0 = -omega * dsdtheta0; 

% 向左偏心(e<0)
% 添加eps避免除零错误
denom = sqrt(L^2 - (r*sin(theta_rad) - e).^2 + eps);
dsdtheta = -r*sin(theta_rad) - ((r*sin(theta_rad)-e).*r.*cos(theta_rad))./denom;
v_left = -omega * dsdtheta; 

% ===== 死点修正 =====
% 在θ=0°和θ=180°强制置零
v0(abs(theta) < 0.5) = 0;
v0(abs(theta-180) < 0.5) = 0;
v_left(abs(theta) < 0.5) = 0;
v_left(abs(theta-180) < 0.5) = 0;

% ===== 绘图 =====
figure;
plot(theta, v0, 'b', 'LineWidth', 1.5);
hold on;
plot(theta, v_left, 'r--', 'LineWidth', 1.5);

% 标注死点位置
plot([0,0], [-200,200], 'k:', 'LineWidth', 1); % TDC线
plot([180,180], [-200,200], 'k:', 'LineWidth', 1); % BDC线
text(-5, 180, 'TDC', 'FontSize', 12, 'HorizontalAlignment', 'right');
text(185, 180, 'BDC', 'FontSize', 12);

% 图表设置
title('滑块速度曲线对比 (含死点修正)');
xlabel('曲柄转角 θ (°)');
ylabel('滑块速度 (mm/s)');
legend('无偏心 (e=0)', '向左偏心 (e=-30mm)', 'Location', 'southeast');
grid on;
xlim([-10, 370]);
xticks(0:90:360);

% 标记关键点
plot(0, 0, 'ko', 'MarkerSize', 10, 'LineWidth', 2); 
plot(180, 0, 'ko', 'MarkerSize', 10, 'LineWidth', 2);
text(75, 120, '工作行程', 'FontSize', 12, 'Color', [0.2 0.5 0.2]);
hold off;