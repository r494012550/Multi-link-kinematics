

%滑块位移主程序
w1=-pi/6;
alpha1=0;
for t=(0:0.01:12);
theta1=w1*t+3*pi/2;
A=L1*cos(theta1)-x0;
B=L1*sin(thetal)+y0;
C=(A^2+B^2+L4^2-L2^2)/(2*L4)
theta4=2*atan((B-(A^2+B^2-C^2)^(1/2))/(A-C))
theta2=atan((B+L4*sin(theta4))/(A+L4*cos(theta4)))
theta3=(theta4)+53*pi/45;
theta6=acos((L1*cos(thetal)+L3*cos(thetal3))/L6)
s=L1+L3+L6+L1*sin(theta1)+L3*sin(theta3)-L6*sin(theta6);
end

%滑块速度主程序
w1=-pi/6;
w2=(L1*w1*((cos(thetal))*(tan(theta4))-sin(thetal)))/(L2*((cos(theta2))*(tan(theta4))-sin(theta2)));
w4=(L1*w1*((cos(theta1))*(tan(theta2))-sin(theta1)))/(L4*(sin(theta4)-(cos(theta4))*(tan(theta2))));
w3=w4;
w6=(L6*sin(theta6))\(L1*wl*sin(thetal)+l3*w3*sin(theta3));
v=L1*wl*cos(thetal)+L3*w3*cos(theta3)-L6*w6*cos(theta6);

%滑块加速度主程序
alpha1=0;
alpha4=(L1*(alphal)*(cos(theta1))*(tan (theta2))-L1*(alphal)*sin(theta1)-L1*w1^2*(sin(thetal))*(tan(theta2))-L1*w1^2*cos(theta1)+L2*w22*sec(theta2)-L4*w4^2*(sin(theta4))*(tan(theta2))-L4*w4^2*cos(theta4))/(L4*sin(theta4)-L4*(cos(theta4))*(tan(theta2)));
alpha3=alpha4;
alpha6=(L1*(alphal)*sin(thetal)+L1*w1^2*cos(thetal)+L3*(alpha3) *sin(theta3)+L3*w3^2*cos(theta3)-L6*w6^2*cos(theta6))/(L6*sin(theta6));a=L1*(alphal)*cos(thetal)-L1*w1^2*sin(theta1)+L3*(alpha3)*cos (theta3)-L3*w3^2*sin(theta3)-L6*(alpha6)*cos (theta6)+L6*w6^2*sin(theta6)