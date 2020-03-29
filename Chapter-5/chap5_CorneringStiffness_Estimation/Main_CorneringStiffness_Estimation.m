function [sys,x0,str,ts] =Main_CorneringStiffness_Estimation(t,x,u,flag)
%***************************************************************%
% Simulation for estimation of cornering stiffness using RLS method
% Assume: lf, lr are known; 
% Assume: Vy, Vel, steering_angle and side-slip angle are measurable
%---------------------------------------------------------------%
% Published by: Kai Liu
% Email:leoking1025@bit.edu.cn
% My github: https://github.com/leoking99-BIT
%***************************************************************% 
    switch flag,
        case 0 % Initialization %
            [sys,x0,str,ts] = mdlInitializeSizes; % Initialization
        case 2 % Update %
            sys = mdlUpdates(t,x,u); % Update discrete states
        case 3 % Outputs %
            sys = mdlOutputs(t,x,u); % Calculate outputs
        case {1,4,9} % Unused flags
            sys = [];            
        otherwise % Unexpected flags %
            error(['unhandled flag = ',num2str(flag)]); % Error handling
    end %  end of switch    
%  End sfuntmpl

function [sys,x0,str,ts] = mdlInitializeSizes
%==============================================================
% Initialization, flag = 0，mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function
%==============================================================
sizes = simsizes;%用于设置模块参数的结构体用simsizes来生成
sizes.NumContStates  = 0;  %模块连续状态变量的个数
sizes.NumDiscStates  = 6;  %模块离散状态变量的个数,实际上没有用到这个数值
sizes.NumOutputs     = 10; %S函数的输出序列的个数
sizes.NumInputs      = 13; %S函数模块输入变量的个数，即CarSim的输出量
sizes.DirFeedthrough = 1;  %模块是否存在直接贯通(direct feedthrough). 
sizes.NumSampleTimes = 1;  %模块的采样次数，>=1
sys = simsizes(sizes);    %设置完后赋给sys输出

x0 = zeros(sizes.NumDiscStates,1);%initial the  state vector
str = [];             % 保留参数，Set str to an empty matrix.
ts  = [0.05 0];       % ts=[period, offset].采样周期sample time=0.05s 

%------------Global parameters and initialization--------------------%
    global InitialGapflag; 
    InitialGapflag = 0; % Ignore the first few inputs from CarSim
    % vehicle parameters initialization
    global VehicleParams; % for SUV
    VehicleParams.Lf  = 1.12;  % 1.05
    VehicleParams.Lr  = 1.48;  % 1.55

    
    % RLS initialization
    [y, c] = func_RLSEstimation_Cf('initial', 0.95, 1, 10);
    [y, c] = func_RLSEstimation_Cr('initial', 0.95, 1, 10);
%     
%     global RLS_params; 
%     RLS_params.nDataTuple    = 10;
%     RLS_params.nCoefficients = 1;
%     RLS_params.initialCoefficients = -90000* ones(RLS_params.nCoefficients,1);    
%     RLS_params.delta         = 10;
%     RLS_params.lambda = 0.95;
%     [y, e, c] = func_RLS_Alt_New('initial', 0, params);   
% 
%     global pre_regressor;
%     nc = 1;
%     nd = 10;
%     pre_regressor.input = zeros(nc,nd);
%     pre_regressor.d     = zeros(nd,1);
    
%  End of mdlInitializeSizes

function sys = mdlUpdates(t,x,u)
%==============================================================
% Update the discrete states, flag = 2， mdlUpdate
% Handle discrete state updates, sample time hits, and major time step
% requirements.
%==============================================================
%  目前没有用到这个过程；在后期的程序模块化时可以继续开发这个功能。
    sys = x;    
% End of mdlUpdate.

function sys = mdlOutputs(t,x,u)
%==============================================================
% Calculate outputs, flag = 3， mdlOutputs
% Return the block outputs. 
%==============================================================
global InitialGapflag;
global VehicleParams;

Cf_Hat        = 0;
Cr_Hat        = 0;
Rinv_f        = 0; 
Rinv_r        = 0;
Hat_err_f     = 0;
Hat_err_r     = 0;  
alpha_f_Hat    = 0;
alpha_r_Hat    = 0; 
alpha_f_Direct = 0;
alpha_r_Direct = 0;

if InitialGapflag < 3 %  get rid of the first two inputs,  because no data from CarSim
    InitialGapflag = InitialGapflag + 1;
else % start control
    InitialGapflag = InitialGapflag + 1;
    %-----Update State Estimation of measured Vehicle Configuration------%
    [Carsim_export] = func_CarsimData_Parse(u);   
    Vx            = Carsim_export.x_dot; 
    Vy            = Carsim_export.y_dot; 
    yawrate       = Carsim_export.phi_dot; % rad/s
    fwa           = Carsim_export.fwa;
    Fyf_Direct    = Carsim_export.Fyf;
    Fyr_Direct    = Carsim_export.Fyr;   
    alpha_f_Direct = Carsim_export.alphaf;
    alpha_r_Direct = Carsim_export.alphar;

    %-----Estimate Cornering stiffness use estimated sideslip angle-----%  
    %for front tire
    alpha_f_Hat = (Vy + yawrate*VehicleParams.Lf)/Vx - fwa;
    [Fyf_hat, Cf_Hat, Rinv_f] = func_RLSEstimation_Cf(alpha_f_Hat, Fyf_Direct);
    Hat_err_f = Fyf_hat - Fyf_Direct;
    %for rear tire 
    alpha_r_Hat = (Vy - yawrate*VehicleParams.Lr)/Vx;
    [Fyr_hat, Cr_Hat, Rinv_r] = func_RLSEstimation_Cr(alpha_r_Hat, Fyr_Direct);
    Hat_err_r = Fyr_hat - Fyr_Direct;
    %-Estimate Cornering stiffness use direct sideslip angle from CarSim-%  
    
    
end % end of if Initialflag < 2 % 
    
sys = [Cf_Hat; Cr_Hat; Rinv_f; Hat_err_f; Rinv_r; Hat_err_r; alpha_f_Hat; alpha_r_Hat; alpha_f_Direct; alpha_r_Direct]; % 

% sys = [t_Elapsed; Ax; Ay_G_SM; fwa; Beta; Vel; Vy; yawrate; Roll; Rollrate; vx_hat; vy_hat; yawrate_hat; roll_hat; rollrate_hat; CafHat; CarHat; C_alpha_f_hat_ay; C_alpha_r_hat_ay]; % 
    
%  sys = [Ctrl_SteerSW; CafHat; CarHat; Fyf; Fyr; alphaf; alphar; Arfa_f; Arfa_r];  

% end  %End of mdlOutputs.

%==============================================================
% sub functions
%==============================================================    

%***************************************************************%
% **** State estimation
%***************************************************************%
function [Sparsed_Carsim_Data] = func_CarsimData_Parse(CarsimData)
%***************************************************************%
% Parse data exported from CarSim, 数据顺序与CarSim的输出一致
%***************************************************************%       
    Sparsed_Carsim_Data.x_dot   = CarsimData(1)/3.6; %Unit:km/h-->m/s，保留1位小数  
    Sparsed_Carsim_Data.y_dot   = CarsimData(2)/3.6; %Unit:km/h-->m/s，保留1位小数   
    Sparsed_Carsim_Data.phi_dot = (round(10*CarsimData(3))/10)*pi/180; %Unit：deg/s-->rad/s，保留1位小数      
    Sparsed_Carsim_Data.fwa     = (round(10*0.5*(CarsimData(4)+ CarsimData(5)))/10)*pi/180; % deg-->rad
    Sparsed_Carsim_Data.alphaf     = (round(10*0.5 * (CarsimData(6)+ CarsimData(8)))/10)*pi/180; % deg-->rad，保留1位小数   
    Sparsed_Carsim_Data.alphar     = (round(10*0.5 * (CarsimData(7)+ CarsimData(9)))/10)*pi/180; % deg-->rad，保留1位小数  
    
    Fy_l1      = round(10*CarsimData(10))/10; %Unit:N，保留1位小数  
    Fy_l2      = round(10*CarsimData(11))/10; %Unit:N，保留1位小数  
    Fy_r1      = round(10*CarsimData(12))/10; %Unit:N，保留1位小数  
    Fy_r2      = round(10*CarsimData(13))/10; %Unit:N，保留1位小数  
    Sparsed_Carsim_Data.Fyf  = Fy_l1 + Fy_r1;
    Sparsed_Carsim_Data.Fyr  = Fy_l2 + Fy_r2;    
% end % end of func_StateEstimation


