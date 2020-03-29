function [vars, status] = MPC_Speed_Controller_CVXGEN(kesi, SpeedProfile, MPCParameters)
% Input:
% Kesi = [Vh, Ah]'
% SpeedProfile: cell(30,1),给定的参考速度。
% MPCParameters： MPC相关参数
% 
% Output:
% status.converged: 1-converge, 0-nonconverge 
% vars.u_0, vars.u{i}: optimized acceleration
% vars.x{i}: predicted states
%
% min J = Q*[Vh(k)-SpeedProfile]^2+R*[ah.des(k)]^2+S*[ah.des(k)-ah.des(k-1)]^2
% s.t. kesi(k+1) = A*kesi(k)*B*u(k)
%---------------------------------------------------------------%
% Published by: Kai Liu
% Email:leoking1025@bit.edu.cn
% My github: https://github.com/leoking99-BIT   
%---------------------------------------------------------------%

%****Step(1): longitudinal vehilce model-one step delay******************%
Ts = MPCParameters.Ts; % 50ms
tao = 0.2; %车辆纵向系统惯性延时参数，0.2 for simulation and 0.5 for real-vehicle test
A =  [1      Ts;
      0      1-Ts/tao];
B =  [0;     Ts/tao]; 

%------Or 忽略惯性延时，此时车辆纵向运动可以表示为一个双积分系统
% A =  [1     Ts;
%       0     1];
% B =  [0;    Ts]; 
    
%****Step(2): MPC formulation;CVXGEN solver *****************************%
    settings.verbose    = 0;       % 0-Silence; 1-display
    settings.max_iters  = 25;    %Limits the total iterations
    settings.resid_tol  = 1e-1;
    settings.eps        = 1e1;
    
    params.x_0      = kesi;
    params.ahm      = kesi(2); % 
    params.vref     = SpeedProfile; % 
    
    params.An      = A;
    params.Bn      = B;
    params.C       = [1, 0];
    
    params.Q        = MPCParameters.Q; %100 
    params.R        = MPCParameters.R; %0.1   
    params.S        = MPCParameters.S; %0.1  

    params.dumax   = MPCParameters.dumax;    
    params.umax    = MPCParameters.umax;
    params.umin    = MPCParameters.umin;    
    
    [vars, status] = csolve_CC(params, settings);
    
end
