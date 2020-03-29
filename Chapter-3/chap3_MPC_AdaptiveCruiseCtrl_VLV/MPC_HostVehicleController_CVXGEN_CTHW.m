function [vars, status] = MPC_HostVehicleController_CVXGEN_CTHW(kesi_host, Ah, Avl_des, Xmin_host, Xmax_host, thw)
% Input:
% kesi_host = [dr_vl; Vr_vl; EgoV.Vh; EgoV.acc ]
% Avl_des £º virtual vehicle optimal control 
% Ah: EgoV.acc
% Xmin_host=[Dsafe;     0; 0;       0];
% Xmax_host=[10000;     0; Vlimits; 0];  
% 
% Self-defined parameters:
% Ts--control period
% dumax--limits of jerk
% umax--the max of acceleration
% umin--the max of deceleration
% 
% Output:
% status.converged: 1-converge, 0-nonconverge 
% vars.u_0, vars.u{i}: optimized control
% vars.x{i}: predicted states
%
%Constant Time-Headway = thw*Vh+r; r:constant
%---------------------------------------------------------------%
% Published by: Kai Liu
% Email:leoking1025@bit.edu.cn
% My github: https://github.com/leoking99-BIT  
%---------------------------------------------------------------%

%****Step(1): longitudinal vehilce model-one step delay******************%
Ts  = 0.05; % 50ms
Ts2 = Ts*Ts;
tao = 0.1; %0.5
A =  [1     Ts  0   -Ts2/2;
      0     1   0   -Ts;
      0     0   1   Ts;
      0     0   0   1-Ts/tao ];
B1 = [0;        0;     0;      Ts/tao]; 
B2 = [Ts2/2;    Ts;    0;      0]; 
C =[1  0   -thw     0];
H =[1 0 0 0;
    0 0 1 0];

% end  
%****Step(2): Constraints and bounds*************************************%
dumax   = 5; % limits of jerk
umax    = 1.5; % the max of acceleration
umin    = -5;% the max of deceleration

%****Step(3): Weighting coeffcients**************************************%
Q = 10;
R = 0.1;
S = 0.5;

%****Step(4): MPC formulation;CVXGEN solver *****************************%
    settings.verbose    = 0;       % 0-Silence; 1-display
    settings.max_iters  = 25;    %Limits the total iterations
    settings.resid_tol  = 1e-1;
    settings.eps        = 1e1;
    
    params.x_0     = kesi_host;
    params.ahm     = Ah; % Ah
    params.avlp    = Avl_des; % Avl_des
    
    params.An      = A;
    params.B1      = B1;
    params.B2      = B2;
    params.C       = C;
    
    params.Q   = Q;
    params.R    = R;    
    params.S    = S;   
    params.r    = 15;%1.5; 
    
    params.dumax   = dumax;    
    params.umax    = umax;
    params.umin    = umin;    

    params.H        = H;   
    params.Xmin     = Xmin_host;
    params.Xmax     = Xmax_host;
    
    [vars, status] = csolve_CTHW(params, settings);
    
end
