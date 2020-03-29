
close all
% firstline          = u.signals.values(:,1);
% Num = length(firstline);

Vx_hat          = u.signals.values(:,1);
Vy_hat          = u.signals.values(:,2);
Yawrate_hat     = u.signals.values(:,3);
Roll_hat        = u.signals.values(:,4);
Rollrate_hat    = u.signals.values(:,5);

Vx_Carsim       = u.signals.values(:,6);
Vy_Carsim       = u.signals.values(:,7);
Yawrate_Carsim  = u.signals.values(:,8);
Roll_Carsim     = u.signals.values(:,9);
Rollrate_Carsim = u.signals.values(:,10);

Ts = 0.05;
Num = length(Vx_hat);
Rollrate_Est = zeros(Num,1);
for i = 2:1:Num
    Rollrate_Est(i) = (Roll_hat(i)-Roll_hat(i-1))/Ts;
end

Vy_error = (Vy_hat - Vy_Carsim); 
Roll_error = (Roll_hat - Roll_Carsim); 
Rollrate_error = (Rollrate_Est - Rollrate_Carsim); 
Rollrate_error1 = (Rollrate_hat - Rollrate_Carsim); 


Estimation_Result.Vy_hat            = Vy_hat;
Estimation_Result.Vy_Carsim         = Vy_Carsim;
Estimation_Result.Roll_hat          = Roll_hat;
Estimation_Result.Roll_Carsim       = Roll_Carsim;
Estimation_Result.Rollrate_Est      = Rollrate_Est;
Estimation_Result.Rollrate_Carsim   = Rollrate_Carsim;
Estimation_Result.Vy_error          = Vy_error;
Estimation_Result.Roll_error        = Roll_error;
Estimation_Result.Rollrate_error    = Rollrate_error;

% save DynamicsEstimation_DualUKF_Modified_10.mat Estimation_Result;

mean_vy         = mean(Vy_error);
var_vy          = var(Vy_error);
mean_Roll       = mean(Roll_error);
var_Roll        = var(Roll_error);
mean_Rollrate   = mean(Rollrate_error);
var_Rollrate    = var(Rollrate_error);
Sumary = [mean_vy, var_vy, mean_Roll, var_Roll, mean_Rollrate, var_Rollrate];

%% figure plot
lw = 2;


figure
plot(1:Num, Vx_hat,'b',1:Num, Vx_Carsim,'r','Linewidth',lw);
legend('Vx_{hat}','Vx_{Carsim}');
title('Comparison of Vx');

figure
plot(1:Num, Yawrate_hat,'b',1:Num, Yawrate_Carsim,'r','Linewidth',lw);
legend('Yawrate_{hat}','Yawrate_{Carsim}');
title('Comparison of Yawrate');

figure
subplot(2,1,1)
plot(1:Num, Roll_hat,'b',1:Num, Roll_Carsim,'r','Linewidth',lw);
legend('Roll_{hat}','Roll_{Carsim}');
title('Comparison of Roll');
subplot(2,1,2)
plot(1:Num, Roll_error,'b','Linewidth',lw);

figure
subplot(2,1,1)
plot(1:Num, Rollrate_hat,'b',1:Num, Rollrate_Carsim,'r','Linewidth',lw);
legend('Rollrate_{hat}','Rollrate_{Carsim}');
title('Comparison of Rollrate');
 subplot(2,1,2)
plot(1:Num, Rollrate_error1,'b',1:Num, Rollrate_error,'k','Linewidth',lw);

figure
subplot(2,1,1)
plot(1:Num, Vy_hat,'b',1:Num, Vy_Carsim,'r','Linewidth',lw);
legend('Vy_{hat}','Vy_{Carsim}');
subplot(2,1,2)
plot(1:Num, Vy_error,'k','Linewidth',lw);
title('Comparison of Vy');




