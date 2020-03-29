
% sys = [Ctrl_SteerSW; t_Elapsed; PosX; PosY; PosPsi; Station; Vel; e_psi; e_d; fwa_opt; fwa_measured]; %
close all
lw = 2;
Num = length(u.signals.values(:,2));

Ctrl_SteerSW     = u.signals.values(:,1);
t_Elapsed   = u.signals.values(:,2);
PosX        = u.signals.values(:,3);
PosY        = u.signals.values(:,4);
PosPsi      = u.signals.values(:,5);
Station     = u.signals.values(:,6);
Vel         = u.signals.values(:,7);
e_psi       = u.signals.values(:,8);
e_y         = u.signals.values(:,9);
fwa_opt          = u.signals.values(:,10);
fwa_measured     = u.signals.values(:,11);

%%
figure % Vel
plot(1:Num, fwa_opt,'b',1:Num, fwa_measured,'k','Linewidth',lw);
grid on
legend('fwa_opt','fwa_measured');

figure % 
plot(1:Num, e_psi,'b','Linewidth',lw);
grid on
legend('e_{psi}');

figure % 
plot(1:Num, e_y,'k','Linewidth',lw);
grid on
legend('e_d');

figure % 
plot(1:Num, t_Elapsed,'k','Linewidth',lw);
grid on
legend('t_Elapsed');


figure % 
plot(PosX, PosY,'b','Linewidth',lw);
grid on
legend('traj');


