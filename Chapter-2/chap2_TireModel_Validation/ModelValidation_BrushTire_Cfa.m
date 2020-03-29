% function ModelValidation_BrushTire_Cfa(u)
%—È÷§Ca = By*Cy*Dy
mu = 1.0;
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
FyL2 = zeros(Num,1);
FyL2C = zeros(Num,1);
FyR2 = zeros(Num,1);
FyR2C = zeros(Num,1);

for index = 1:1:Num
    [Calphar1(index),Car1(index),FyL2(index),FyL2C(index)] = func_BrushTyreModel_puresideslip(Fz_L2(index), alpha_L2(index), mu);
    [Calphar2(index),Car2(index),FyR2(index),FyR2C(index)] = func_BrushTyreModel_puresideslip(Fz_R2(index), alpha_R2(index), mu);    
    CarM1(index) = Fy_L2(index)/(alpha_L2(index)*pi/180);
    CarM2(index) = Fy_R2(index)/(alpha_R2(index)*pi/180);
    
    FyL2MF(index) = func_MFTyreModel_puresideslip(Fz_L2(index), alpha_L2(index));
    FyR2MF(index) = func_MFTyreModel_puresideslip(Fz_R2(index), alpha_R2(index));
    
end

Fyr   = Fy_L2 + Fy_R2;
FyrC   = FyL2C + FyR2C;
FyrBT  = FyL2 + FyR2;
FyrMF  = FyL2MF + FyR2MF;

figure (1)
plot(iStart:Num, Fyr,'B',iStart:Num, FyrBT,'b*',iStart:Num, FyrC,'r',iStart:Num, FyrMF,'r*');
grid on

% figure (1)
% plot(iStart:Num, Fy_L2,'B',iStart:Num, FyL2,'b*',iStart:Num, FyL2C,'r');
% grid on
% figure (2)
% plot(iStart:Num, Fy_R2,'B',iStart:Num, FyR2,'b*',iStart:Num, FyR2C,'r');
% grid on

% figure (3)
% plot(iStart:Num, Calphar1,'b',iStart:Num, Car1,'b*',iStart:Num, CarM1,'r');
% grid on
% figure (4)
% plot(iStart:Num, Calphar2,'b',iStart:Num, Car2,'b*',iStart:Num, CarM2,'r');
% grid on

% end % End of func