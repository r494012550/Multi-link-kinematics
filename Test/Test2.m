% 参数
R = 0.5;    % 曲柄长度
L = 2;      % 连杆长度
theta = linspace(0, 2*pi, 100); % 曲柄转角

% 复数法计算
phi = -asin(R/L * sin(theta));  % 从虚部方程解出phi
x = R*cos(theta) + L*cos(phi);  % 代入实部方程

% 可视化
figure;
plot(theta, x, 'LineWidth', 2);
xlabel('曲柄转角 \theta (rad)');
ylabel('滑块位移 x (m)');
title('复数法计算的滑块位移曲线');
grid on;