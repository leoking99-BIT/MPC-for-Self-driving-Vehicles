%***************************************************************%
% LinearTireModel_WithKs_FOH
% Deem Vx as constant
%***************************************************************%
function [StateSpaceModel] = func_RigidbodyDynamicalModel_FOH_Extended(VehiclePara, MPCParameters, Vel, Calpha_f, Calpha_r)
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

K_phi       = 145330;
D_phi       = 4500;
Theta_1     = Calpha_f + Calpha_r;
Theta_2     = lf*Calpha_f - lr*Calpha_r;
Theta_3     = lf*lf*Calpha_f + lr*lr*Calpha_r;

Gai1 = hCG*Theta_2/(Ix*Vel)- hCG*M*Vel/Ix;
Acn = [ Theta_1/(M*Vel)         Theta_2/(M*Vel)-Vel     0                0                       0   0;
        Theta_2/(Iz*Vel)        Theta_3/(Iz*Vel)        0                0                       0   0;
        hCG*Theta_1/(Ix*Vel)    Gai1                    -D_phi/Ix       (hCG*M*g-K_phi)/Ix     0   0;
        0                       0                       1               0                       0   0;
        1                       0                       0               0                       0   Vel;
        0                       1                       0               0                       0   0 ];
    
B1cn = [-Calpha_f/M;	-lf*Calpha_f/Iz;    -hCG*Calpha_f/Ix;    0;     0; 	0]; % 
B2cn = [-g,              0,                 -M*g*hCG/Ix,         0,     0, 	0;
        0,               0,                  0,                  0,     0, 	-Vel]';

% Acn = [ Theta_1/(M*Vel)         Theta_2/(M*Vel)-Vel     0                0                       0   0;
%         Theta_2/(Iz*Vel)        Theta_3/(Iz*Vel)        0                0                       0   0;
%        0                        0                       0               0                       0    0;
%         0                       0                       0               0                       0   0;
%         1                       0                       0               0                       0   Vel;
%         0                       1                       0               0                       0   0 ];
%     
% B1cn = [-Calpha_f/M;	-lf*Calpha_f/Iz;    0;          0;     0; 	0]; % 
% B2cn = [-g,              0,                 0,         0,     0, 	0;
%         0,               0,                  0,        0,     0, 	-Vel]';
    
    As = eye(Nx) + Ts * Acn;
    Bs1 = Ts * B1cn; 
    Bs2 = Ts * B2cn;
    
    AcnTsl = Tsl*Acn;
    Al = eye(Nx) + AcnTsl + 0.5*AcnTsl*AcnTsl ; %  + (AcnTsl*AcnTsl*AcnTsl)/6;

    Baug    = [B1cn, B2cn];
    BaugTsl = Baug*Tsl;
    Gail1   = BaugTsl + 0.5*(AcnTsl * BaugTsl); % + (AcnTsl * AcnTsl * BaugTsl)/6;
    Gail2   = 0.5*BaugTsl;% + (AcnTsl * BaugTsl)/6;
    
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

