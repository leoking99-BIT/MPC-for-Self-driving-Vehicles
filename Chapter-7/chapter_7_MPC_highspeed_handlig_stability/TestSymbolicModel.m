
syms K_phi D_phi Theta_1 Theta_2 Theta_3
syms M hCG Iz Ix Vel g Calpha_f Calpha_r lf


Mint = [  M         0       -M*hCG          0;
          0         Iz      0               0;
          -M*hCG    0       Ix+M*hCG*hCG    0
          0         0       0               1];

Nint = [-Theta_1/Vel    M*Vel-Theta_2/Vel   0           M*g;  %0-->M*g
        -Theta_2/Vel    -Theta_3/Vel        0           0;
        0               -M*hCG*Vel          D_phi       K_phi-M*g*hCG;
        0               0                   -1          0];

F1int = [-Calpha_f; -lf*Calpha_f;  0;  0];
F2int = [   0               0; 
            0               0;
            K_phi           0;
            0               0];

Ac_11     = -Mint\Nint; % 4*4
B1cn_11   = Mint\F1int; % 4*1
B2cn_11   = Mint\F2int; % 4*2

% K_phi       = 145330;
% D_phi       = 4500;
% Theta_1     = Calpha_f + Calpha_r;
% Theta_2     = lf*Calpha_f - lr*Calpha_r;
% Theta_3     = lf*lf*Calpha_f + lr*lr*Calpha_r;
% 
% Gai1 = hCG*Theta_2/(Ix*Vel)- hCG*M*Vel/Ix;
% Acn = [ Theta_1/(M*Vel)         Theta_2/(M*Vel)-Vel     0                0                       0   0;
%         Theta_2/(Iz*Vel)        Theta_3/(Iz*Vel)        0                0                       0   0;
%         hCG*Theta_1/(Ix*Vel)    Gai1                    -D_phi/Ix       (hCG*M*g-K_phi)/Ix     0   0;
%         0                       0                       1               0                       0   0;
%         1                       0                       0               0                       0   Vel;
%         0                       1                       0               0                       0   0 ];
%     
% B1cn = [-Calpha_f/M;	-lf*Calpha_f/Iz;    -hCG*Calpha_f/Ix;    0;     0; 	0]; % 
% B2cn = [-g,              0,                 -M*g*hCG/Ix,         0,     0, 	0;
%         0,               0,                  0,                  0,     0, 	-Vel]';







