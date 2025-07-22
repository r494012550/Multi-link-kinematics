% 测试两种机构机械增益
L_toggle = [0.15, 0.65, 0.35, 1.1]; % 优化后肘杆参数
L_crank = [0.3, 0.7]; % 曲柄连杆参数

theta = linspace(0, 2*pi, 100);
MA_toggle = arrayfun(@(x) calculate_MA(x, L_toggle), theta);
MA_crank = 0.7*sin(theta); % 曲柄连杆理论增益

plot(theta, MA_toggle, 'r-', theta, MA_crank, 'b--');
legend('肘杆机构','曲柄连杆');
xlabel('驱动角(rad)'); ylabel('机械增益');
title('机械增益对比');
grid on;

% 肘杆机构向量分析法核心代码
function [theta2, theta3, Fout] = toggle_mechanism(theta1, L, Fin)
    % 参数：L = [L1, L2, L3, L4]
    % 向量分解
    AB = L(1)*[cos(theta1); sin(theta1)];
    BC = L(2)*[cos(theta2); sin(theta2)];
    CD = L(3)*[cos(theta3); sin(theta3)];
    AD = [L(4); 0];
    
    % 闭合方程：AB + BC + CD = AD
    eq1 = AB(1) + BC(1) + CD(1) - AD(1);
    eq2 = AB(2) + BC(2) + CD(2);
    
    % 牛顿迭代求解角度
    options = optimset('Display','off');
    sol = fsolve(@(x) [eq1; eq2], [0; pi/2], options);
    theta2 = sol(1);
    theta3 = sol(2);
    
    % 力学分析（虚功原理）
    J = [-L(2)*sin(theta2), -L(3)*sin(theta3);
          L(2)*cos(theta2),  L(3)*cos(theta3)];
    Fout = -J' \ [Fin; 0];
end

% 优化目标：最大化机械增益
function [opt_params, max_gain] = optimize_toggle_mechanism()
    % 初始参数：[L1, L2, L3, L4]
    L0 = [0.2, 0.5, 0.5, 1.0];
    lb = [0.1, 0.3, 0.3, 0.8];
    ub = [0.3, 0.7, 0.7, 1.2];
    
    % 遗传算法优化
    options = optimoptions('ga', 'PopulationSize', 50, 'MaxGenerations', 100);
    [opt_params, max_gain] = ga(@(L) -mechanical_advantage(L), 4, [], [], [], [], lb, ub, [], options);
    max_gain = -max_gain;
    
    function MA = mechanical_advantage(L)
        theta1_range = linspace(pi/6, 5*pi/6, 50);
        MA_values = arrayfun(@(x) calculate_MA(x, L), theta1_range);
        MA = min(MA_values); % 取最小增益保障稳定性
    end

    function ma = calculate_MA(theta1, L)
        [~, ~, Fout] = toggle_mechanism(theta1, L, 1);
        ma = abs(Fout(1));
    end
end

function animate_toggle_mechanism(L, theta1_range)
    figure('Position', [100 100 800 600]);
    axis equal; grid on; hold on;
    xlim([-0.5, L(4)+0.5]);
    ylim([-0.5*L(4), 0.5*L(4)]);
    
    % 初始化图形对象
    h_AB = line([0,0], [0,0], 'Color','b', 'LineWidth',2);
    h_BC = line([0,0], [0,0], 'Color','r', 'LineWidth',2);
    h_CD = line([0,0], [0,0], 'Color','g', 'LineWidth',2);
    h_D = plot(0,0, 'ko', 'MarkerFaceColor','k', 'MarkerSize',8);
    
    % 动画循环
    for theta1 = theta1_range
        [theta2, theta3] = toggle_mechanism(theta1, L, 0);
        
        % 计算关键点坐标
        A = [0; 0];
        B = L(1)*[cos(theta1); sin(theta1)];
        C = B + L(2)*[cos(theta2); sin(theta2)];
        D = [L(4); 0];
        
        % 更新图形数据
        set(h_AB, 'XData', [A(1), B(1)], 'YData', [A(2), B(2)]);
        set(h_BC, 'XData', [B(1), C(1)], 'YData', [B(2), C(2)]);
        set(h_CD, 'XData', [C(1), D(1)], 'YData', [C(2), D(2)]);
        set(h_D, 'XData', D(1), 'YData', D(2));
        
        drawnow;
        pause(0.05);
    end
end

