
close all
Num = length(u.signals.values(:,1));
Car_X = u.signals.values(:,1); 
Car_Y = u.signals.values(:,2); 
Yaw = u.signals.values(:,3); 
X_pred = u.signals.values(:,4); 
Y_pred = u.signals.values(:,5); 
Yaw_pred = u.signals.values(:,6); 
Vx_m_s = u.signals.values(:,7); 
Steer_rad = u.signals.values(:,8); 
Steer_SW = u.signals.values(:,9); 
Steer_ratio = u.signals.values(:,10);


ratio = mean(Steer_ratio);
fprintf('the mean-ratio between Steer_SW and front-wheel-steering is: %f\n', ratio); 

%%
lw=2;
figure
plot(Car_X, Car_Y,'b',X_pred, Y_pred,'r','Linewidth',lw);
legend('CarSim Pos','ModelPredict Pos');
xlabel('X (m)','FontName','Times New Roman','FontSize',14)
ylabel('Y(m)','FontName','Times New Roman','FontSize',14,'Rotation',90)
title('Pos compare');
axis tight

figure
plot(1:Num, Yaw,'b',1:Num, Yaw_pred,'r','Linewidth',lw);
legend('Yaw','Yaw pred');
xlabel('Data Index','FontName','Times New Roman','FontSize',14)
ylabel('Yaw angle (rad)','FontName','Times New Roman','FontSize',14,'Rotation',90)
title('Yaw angle predict');
axis tight

figure
plot(1:Num, Steer_SW, 'b',1:Num, (ratio*180/pi)*Steer_rad,'k', 'Linewidth',lw);
grid on
legend('Steer SW','Desired steer SW');
xlabel('Data Index','FontName','Times New Roman','FontSize',14)
ylabel('Steer_SW angle (deg)','FontName','Times New Roman','FontSize',14,'Rotation',90)
title('Steer SW predict');
axis tight

