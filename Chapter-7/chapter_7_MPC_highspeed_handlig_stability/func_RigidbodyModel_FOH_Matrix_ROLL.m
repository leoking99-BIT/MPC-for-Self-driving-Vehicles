%***************************************************************%
% LinearTireModel_WithKs_FOH
% Deem Vx as constant
%***************************************************************%
function [StateSpaceModel] = func_RigidbodyModel_FOH_Matrix_ROLL(VehiclePara, MPCParameters, Vel, Calpha_f, Calpha_r)
    % generate State-space model
   
M       = VehiclePara.m ;
g       = VehiclePara.g;%m/s^2
L       = VehiclePara.L;  %a = 1.11;  
lf      = VehiclePara.Lr;% 1.12;%1.05%  重心向后转移了0.7m
lr      = VehiclePara.Lr;% 1.48;%
Tr      = VehiclePara.Tr ; %m
hCG     = VehiclePara.hCG;%m
Ix      = VehiclePara.Ix;
Iz      = VehiclePara.Iz;

Np = MPCParameters.Np;
Ns = MPCParameters.Ns;
Nx = MPCParameters.Nx;
Ny = MPCParameters.Ny;
Ts = MPCParameters.Ts;
Tsl = MPCParameters.Tsl; 

% Np      = 25; % 0.1
% Ns      = 10; % 0.05;
% Nx      = 8;
% Ny      = 5;
% Ts      = 0.05;
% Tsl     = 0.5;

%% 
K_phi       = 145330;
D_phi       = 4500;
Theta_1     = Calpha_f + Calpha_r;
Theta_2     = lf*Calpha_f - lr*Calpha_r;
Theta_3     = lf*lf*Calpha_f + lr*lr*Calpha_r;

Ixz     = 0; %500
Mint = [  M         0       -M*hCG          0;
          0         Iz      -Ixz            0;
          -M*hCG    -Ixz    Ix+M*hCG*hCG    0
          0         0       0               1];

Nint = [-Theta_1/Vel    M*Vel-Theta_2/Vel   0           M*g;  %0-->M*g
        -Theta_2/Vel    -Theta_3/Vel        0           0;
        0               -M*hCG*Vel          D_phi       K_phi-M*g*hCG;
        0               0                   -1          0];

F1int = [-Calpha_f; -lf*Calpha_f;  0;  0];
F2int = [   0        0; 
            0        0;
            K_phi    0; %-M*g*hCG
            0        0];
    
% M_inv     = inv(M);
% Ac_11     = M_inv*N;
Ac_11     = -Mint\Nint; % 4*4
B1cn_11   = Mint\F1int; % 4*1
B2cn_11   = Mint\F2int; % 4*2

Ac_12     = [0 0; 0 0; 0 0; 0 0];
Ac_21     = [1  0    0    0;
             0  1    0    0];
Ac_22     = [0   Vel; 0   0]; 
 
B1cn_12   = [ 0; 	0]; 
B2cn_12   = [ 0, 	0; 
              0, 	-Vel];
    
Acn     = [Ac_11 Ac_12; Ac_21   Ac_22];
B1cn    = [B1cn_11; B1cn_12];
B2cn    = [B2cn_11; B2cn_12];

%% ZOH/FOH
    As = eye(Nx) + Ts * Acn;
    Bs1 = Ts * B1cn; 
    Bs2 = Ts * B2cn;
    
    AcnTsl = Tsl*Acn;
    Al = eye(Nx) + AcnTsl + 0.5*AcnTsl*AcnTsl  + (AcnTsl*AcnTsl*AcnTsl)/6; %; % 

    Baug    = [B1cn, B2cn];
    BaugTsl = Baug*Tsl;
    Gail1   = BaugTsl + 0.5*(AcnTsl * BaugTsl)+ (AcnTsl * AcnTsl * BaugTsl)/6;%; % 
    Gail2   = 0.5*BaugTsl + (AcnTsl * BaugTsl)/6;%;%
    
    B1l     = Gail1 - Gail2; 
    B2l     = Gail2;
    Bl11    = B1l(:,1);
    Bl12    = B1l(:,2:3);
    Bl21    = B2l(:,1);
    Bl22    = B2l(:,2:3);
    
    StateSpaceModel.Acn     = Acn;
    StateSpaceModel.B1cn    = B1cn; 
    StateSpaceModel.B2cn    = B2cn; 
    
    StateSpaceModel.As      = As;
    StateSpaceModel.Bs1     = Bs1; 
    StateSpaceModel.Bs2     = Bs2;
    
    StateSpaceModel.Al      = Al;
    StateSpaceModel.Bl11    = Bl11; 
    StateSpaceModel.Bl12    = Bl12; 
    StateSpaceModel.Bl21    = Bl21; 
    StateSpaceModel.Bl22    = Bl22; 
    
end % end of func_SpatialDynamicalModel

