% a = 1 + 2i;
% b = 3 + 4i;
% c = a + b;
% d = a / b;
%  disp(c);
%  disp(d);

% a = 2;
% b = 2;
% c = 2*sqrt(2);
% syms theta1 real;
% syms theta2 real;
% syms theta3 real;

% eq = c*exp(1i*theta3) == a*exp(1i*theta1) + b*exp(1i*theta2);

% eq_r = real(eq);
% eq_i = imag(eq);

% % 求解方程组
% solution = solve([eq_r, eq_i ,  theta1 == pi/4], [theta1, theta2, theta3], 'Real', true);

% disp('theta2 = '); disp(double(solution.theta2))
% disp('theta3 = '); disp(double(solution.theta3))

% a = 2; b = 2; c = 2*sqrt(2);
% syms theta1 theta2 theta3 real

% % 复数向量方程
% eq = c*exp(1i*theta3) == a*exp(1i*theta1) + b*exp(1i*theta2);

% % 分离实部虚部
% eq_real = real(eq) == 0;
% eq_imag = imag(eq) == 0;

% % 添加驱动约束：theta1为已知输入（例如pi/4）
% theta1_val = pi/4;  
% constraint = [eq_real, eq_imag, theta1 == theta1_val];

% % 求解剩余两个角度
% solution = solve(constraint, [theta2, theta3], 'Real', true);
% disp('theta2 = '); disp(double(solution.theta2))
% disp('theta3 = '); disp(double(solution.theta3))

a = 2; b = 2; c = 2*sqrt(2);
theta1_val = linspace(0, 2*pi, 100); % 曲柄旋转一周

for i = 1:length(theta1_val)
    % 解非线性方程组
    fun = @(x) [
        a*cos(theta1_val(i)) + b*cos(x(1)) - c*cos(x(2));
        a*sin(theta1_val(i)) + b*sin(x(1)) - c*sin(x(2))
    ];
    
    % 初始值[0,0], 求解theta2和theta3
    x0 = [0, 0];
    options = optimset('Display','off');
    sol = fsolve(fun, x0, options);
    
    theta2(i) = sol(1);
    theta3(i) = sol(2);
end

% 可视化
plot(theta1_val, theta2, 'r', theta1_val, theta3, 'b');
legend('\theta_2 (连杆角)', '\theta_3 (输出角)');
xlabel('\theta_1 (曲柄角)');
grid on