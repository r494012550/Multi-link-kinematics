% 机构参数初始化
params.L = [237, 1150,1020, 800, 1000]; % [L1, L2, L3, L4, L6]
params.origin = [0; 0];        % [x0; y0]
w1 = -pi/6;                    % 曲柄角速度
t = 0:0.01:12;                 % 时间序列

% 向量化计算核心（修正所有希腊字母变量名）
theta1 = w1*t + 3*pi/2;        % 曲柄角向量

% 复数向量表示各关键点（修正参数命名）
P1 = params.L(1)*exp(1i*theta1); 
D = params.L(4)*exp(1i*asin((imag(P1)+params.origin(2))/params.L(4))); 

% 闭合矢量方程求解（统一使用theta命名）
P2 = real(P1) + 1i*(imag(P1) + params.origin(2)) - D;
theta2 = angle(P2 - params.origin(1)); 
theta3 = angle(D) + 53*pi/45; 

% 滑块位置计算（修正lambda函数参数名）
slider_pos = @(theta6) params.L(5)*exp(1i*theta6);
eq = @(theta6) params.L(1)*cos(theta1) + params.L(3)*cos(theta3) - real(slider_pos(theta6));
theta6 = acos(eq(0)/params.L(5));  % 统一使用theta6

% 运动计算结果
s = params.L(1)+params.L(3)+params.L(5) + imag(P1) + ...
    params.L(3)*sin(theta3) - params.L(5)*sin(theta6);

% 可视化
figure;
subplot(3,1,1), plot(t,s), title('位移'), xlabel('t/s'), ylabel('s/m')
subplot(3,1,2), plot(t,v), title('速度'), xlabel('t/s'), ylabel('v/(m/s)')
subplot(3,1,3), plot(t,a), title('加速度'), xlabel('t/s'), ylabel('a/(m/s²)')