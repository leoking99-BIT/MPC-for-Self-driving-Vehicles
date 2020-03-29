% function ModelValidation_MFTire(u)
%缺点是不能体现mu在模型中的作用。
%跟定的参数是在mu=1的情况下测试的。
clc
close all

mu = 1.0; % 0.3

alpha_L1 = u.signals.values(:,14);
alpha_L2 = u.signals.values(:,15);
alpha_R1 = u.signals.values(:,16);
alpha_R2 = u.signals.values(:,17);
Fz_L1    = u.signals.values(:,18);
Fz_L2    = u.signals.values(:,19);
Fz_R1    = u.signals.values(:,20);
Fz_R2    = u.signals.values(:,21);
Fy_L1    = u.signals.values(:,22);
Fy_L2    = u.signals.values(:,23);
Fy_R1    = u.signals.values(:,24);
Fy_R2    = u.signals.values(:,25);

iStart = 1;
Num = length(alpha_L1);
FyL1 = zeros(Num,1);
FyL2 = zeros(Num,1);
FyR1 = zeros(Num,1);
FyR2 = zeros(Num,1);
for index = 1:1:Num
    FyL1(index) = func_MFTyreModel_puresideslip(Fz_L1(index), alpha_L1(index));
    FyL2(index) = func_MFTyreModel_puresideslip(Fz_L2(index), alpha_L2(index));
    FyR1(index) = func_MFTyreModel_puresideslip(Fz_R1(index), alpha_R1(index));
    FyR2(index) = func_MFTyreModel_puresideslip(Fz_R2(index), alpha_R2(index));
end

figure (1)
plot(iStart:Num, Fy_L1,'b',iStart:Num, FyL1,'b*',iStart:Num, Fy_R1,'r',iStart:Num, FyR1,'r*');
legend('Fyf-L1-CarSim','Fyf-L1-Pacejka','Fyf-R1-CarSim','Fyf-R1-Pacejka');
grid on

figure (2)
plot(iStart:Num, Fy_L2,'b',iStart:Num, FyL2,'b*',iStart:Num, Fy_R2,'r',iStart:Num, FyR2,'r*');
legend('Fyr-L2-CarSim','Fyf-L2-Pacejka','Fyf-R2-CarSim','Fyf-R2-Pacejka');
grid on

% end % End of func