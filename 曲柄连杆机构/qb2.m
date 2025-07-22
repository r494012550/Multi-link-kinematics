% 机构参数
L1 = 160;
L2 = 1230;
theta = 0:1:360;
thetarad = deg2rad(theta);

x = @(th) L1*sin(th) + sqrt(L2^2-(L1*cos(th)).^2);

v = gradient(x(thetarad),thetarad);

a = gradient(v,thetarad);



% 分别绘制位置、速度和加速度
figure;
subplot(3,1,1);
plot(theta, x(thetarad));
title('滑块位置随时间变化');
xlabel('角度(θ)');
ylabel('位置 (mm)');

subplot(3,1,2);
plot(theta, v);
title('滑块速度随时间变化');
xlabel('角度(θ)');
ylabel('速度 (mm/θ)');

subplot(3,1,3); 
plot(theta, a);
title('滑块加速度随时间变化');
xlabel('角度(θ)');
ylabel('加速度 (mm/θ²)');