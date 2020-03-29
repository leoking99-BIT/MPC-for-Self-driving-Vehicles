
close all
lw = 2;
Num = length(u.signals.values(:,2));

% LTR = u.signals.values(:,9);
% Yzmp_SH = u.signals.values(:,10);
% Yzmp = u.signals.values(:,11);
% Yzmp1 = u.signals.values(:,12);
% 
% figure
% plot(1:Num, LTR','r',1:Num, Yzmp_SH,'b',1:Num, Yzmp,'g',1:Num, Yzmp1,'k', 'Linewidth',lw);
% grid on
% legend('LTR','Yzmp_{SH}','Yzmp','Yzmp1');

%%

Yzmp_SH = u.signals.values(:,10);
Yzmp = u.signals.values(:,11);
LTR = u.signals.values(:,12);

figure
plot(1:Num, LTR,'r',1:Num, Yzmp_SH,'b',1:Num, Yzmp,'g', 'Linewidth',lw);
grid on
legend('LTR','Yzmp_{SH}','Yzmp');

% Roll_BankR = u.signals.values(:,17);
% Roll_shad  = u.signals.values(:,18);
% figure
% plot(1:Num, Roll_shad,'r',1:Num, Roll_BankR,'bo', 'Linewidth',lw);
% grid on
% legend('Roll_{shad}','Roll_{BankR}');

%% Save simulation data
% SimResult_Both = u.signals.values; 
% save SimResult_Both.mat SimResult_Both

% SimResult_Neither = u.signals.values; 
% save SimResult_Neither.mat SimResult_Neither
% 
% SimResult_Bank = u.signals.values; 
% save SimResult_Bank.mat SimResult_Bank
% 
% SimResult_Curvature = u.signals.values; 
% save SimResult_Curvature.mat SimResult_Curvature

