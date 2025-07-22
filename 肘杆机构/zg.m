% 初始猜测值 (单位: 弧度)
x0 = [0.5, -0.5, 0, 0, pi, 1000]; % 示例初始值: [θ2≈28.6°, θ3≈-28.6°, θ4=0°, θ5=0°, θ6=180°, s=1000]

% 设置求解选项
options = optimoptions('fsolve', 'Display', 'iter', 'MaxIterations', 1000);

% 调用fsolve求解
[x_sol, fval] = fsolve(@equations, x0, options);

% 将弧度转换为角度 (前5个变量)
theta_deg = rad2deg(x_sol(1:5));
s_value = x_sol(6);

% 输出结果
disp('求解结果 (角度单位为度):');
disp(['θ2 = ', num2str(theta_deg(1))]);
disp(['θ3 = ', num2str(theta_deg(2))]);
disp(['θ4 = ', num2str(theta_deg(3))]);
disp(['θ5 = ', num2str(theta_deg(4))]);
disp(['θ6 = ', num2str(theta_deg(5))]);
disp(['s  = ', num2str(s_value)]);

function F = equations(x)
    % 输入变量: x = [θ2, θ3, θ4, θ5, θ6, s] (单位: 弧度)
    theta2 = x(1); theta3 = x(2); theta4 = x(3); theta5 = x(4); theta6 = x(5); s = x(6);
    
    % 杆长参数
    L1 = 280; 
    L2 = 1588.5; 
    L3 = 940; 
    L4 = 840; 
    L5 = 850; 
    L6 = 845;
    X0 = 2130; 
    Y0 = 350;
    
    % 方程组的6个分量
    F = zeros(6,1);
    F(1) = L1 + L2*cos(theta2) + L3*cos(theta3) - X0; % 实部(方程1)
    F(2) = L2*sin(theta2) + L3*sin(theta3) + Y0;     % 虚部(方程1)
    F(3) = L4*cos(theta4) + L5*cos(theta5) - L3*cos(theta3); % 实部(方程2)
    F(4) = L4*sin(theta4) + L5*sin(theta5) - L3*sin(theta3); % 虚部(方程2)
    F(5) = L4*cos(theta4) + L6*cos(theta6);           % 实部(方程3)
    F(6) = L4*sin(theta4) + L6*sin(theta6) + s;       % 虚部(方程3)
end