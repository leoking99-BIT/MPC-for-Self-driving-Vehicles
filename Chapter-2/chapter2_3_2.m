%---------------------------------------------------------------%
% Published by: Kai Liu
% Email:leoking1025@bit.edu.cn
% My github: https://github.com/leoking99-BIT
%---------------------------------------------------------------%
%%
% 声明模型推导过程中需呀用到的参数变量
syms K_phi D_phi Theta_1 Theta_2 Theta_3
syms M hCG Iz Ix Vel g Calpha_f Calpha_r lf

%根据实际需要定义车辆动力学微分方程：Mint，Nint，F1int，F2int
Mint = [  M         0       -M*hCG          0;
          0         Iz      0               0;
          -M*hCG    0       Ix+M*hCG*hCG    0;
          0         0       0               1];
Nint = [-Theta_1/Vel    M*Vel-Theta_2/Vel   0           M*g;
        -Theta_2/Vel    -Theta_3/Vel        0           0;
        0               -M*hCG*Vel          D_phi       K_phi-M*g*hCG;
        0               0                   -1          0];
F1int = [-Calpha_f; -lf*Calpha_f;  0;  0];
F2int = [   0               0; 
            0               0;
            K_phi           0;
            0               0];
%通过矩阵变换得到状态空间方程形式的车辆动力学方程
% dot{kesi} = Ac_11*kesi + B1cn_11 * u1 + B2cn_11*u2
Ac_11     = -Mint\Nint;
B1cn_11   = Mint\F1int;
B2cn_11   = Mint\F2int;







