% 曲柄滑块机构运动分析
% 参数定义
r = 160;        % 曲柄半径 (mm)
l = 1230;       % 连杆长度 (mm)
omega = 2*pi;   % 曲柄角速度 (rad/s)
dt = 0.001;     % 时间步长 (s)

% 时间向量 (0 到 1 秒，完整周期)
t = 0:dt:1;     % 1秒对应完整周期(ω=2π rad/s)

% 曲柄角度 (随时间变化)
theta = omega * t + pi / 2;

% 滑块位置函数
x = @(th) r*cos(th) + sqrt(l^2 - (r*sin(th)).^2);

% % 数值微分计算速度和加速度
% v = zeros(size(t));     % 滑块速度
% a = zeros(size(t));     % 滑块加速度

% % 使用中心差分法
% for i = 2:length(t)-1
%     % 速度 = dx/dt
%     v(i) = (x(theta(i+1)) - x(theta(i-1))) / (2*dt);
    
%     % 加速度 = dv/dt
%     a(i) = (v(i+1) - v(i-1)) / (2*dt);
% end

%简单求导
v = gradient(x(theta), dt);
a = gradient(v, dt); 

% % 边界点处理 (前向/后向差分)
% v(1) = (x(theta(2)) - x(theta(1))) / dt;
% v(end) = (x(theta(end)) - x(theta(end-1))) / dt;

% a(1) = (v(2) - v(1)) / dt;
% a(end) = (v(end) - v(end-1)) / dt;

% 显示关键时间点的运动参数
fprintf('时间(s)\t角度(rad)\t位置(mm)\t速度(mm/s)\t加速度(mm/s²)\n');
fprintf('-----------------------------------------------------------------\n');
for ti = [0,0.212,0.25,0.5]   %输入时间(s)得到其他参数
    index = round(ti/dt) + 1;
    fprintf('%.3f\t%.2f\t\t%.2f\t\t%.2f\t\t%.2f\n', ...
            t(index), theta(index)/pi*180, x(theta(index)), v(index), a(index));
end