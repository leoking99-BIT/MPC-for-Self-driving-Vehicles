function [sys,x0,str,ts] =Main_CorneringStiffness_Comparison_RLS_Alt_new(t,x,u,flag)
%***************************************************************%
% lf, lr are known;
% ay, ax, delta_f and Beta are measurable
% Using RLS to online identification [C_alpha_f_hat, C_alpha_r_hat]
% 
%---------------------------------------------------------------%
% Published by: Kai Liu
% Email:leoking1025@gmail.com 
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
% Return the sizes, initial conditions, and sample times for the S-function.
%==============================================================
sizes = simsizes;%用于设置模块参数的结构体用simsizes来生成
sizes.NumContStates  = 0;  %模块连续状态变量的个数
sizes.NumDiscStates  = 6;  %模块离散状态变量的个数,实际上没有用到这个数值，只是用这个来表示离散模块
sizes.NumOutputs     = 4;  %S函数的输出，包括控制量和其它监测量
sizes.NumInputs      = 13; %S函数模块输入变量的个数，即CarSim的输出量
sizes.DirFeedthrough = 1;  %模块是否存在直接贯通(direct feedthrough). 1 means there is direct feedthrough.
% 直接馈通表示系统的输出或可变采样时间是否受到输入的控制。
% a.  输出函数（mdlOutputs或flag==3）是输入u的函数。即，如果输入u在mdlOutputs中被访问，则存在直接馈通。
% b.  对于一个变步长S-Function的“下一个采样时间”函数（mdlGetTimeOfNextVarHit或flag==4）中可以访问输入u。
% 正确设置直接馈通标志是十分重要的，因为它影响模型中块的执行顺序，并可用检测代数环。
sizes.NumSampleTimes = 1;  %模块的采样次数，>=1

sys = simsizes(sizes);    %设置完后赋给sys输出

x0 = zeros(sizes.NumDiscStates,1);%initial the  state vector， of no use

str = [];             % 保留参数，Set str to an empty matrix.

ts  = [0.05 0];       % ts=[period, offset].采样周期sample time=0.05,50ms 

%------------Global parameters and initialization--------------------%
    % online identification RLS initialization
    [y, e] = func_RLSFilter_Calpha_f_useFyf('initial', 0.95, 10, 10);
    [y, e] = func_RLSFilter_Calpha_r_useFyr('initial', 0.95, 10, 10);
 
    global params; 
    params.nCoefficients = 1;
    params.delta = 10;
    params.initialCoefficients = [-110000];
    params.lambda = 0.95;
%     params.nDataTuple = 2; %filter order, number of data tuple
    [y, e, c] = func_RLS_Alt_New('initial', 0, params);   

    global pre_regressor;
    nc = 1;
    nd = 10;
    pre_regressor.input = zeros(nc,nd);
    pre_regressor.d     = zeros(nd,1);
    
    % vehicle parameters initialization
    global InitialGapflag; 
    InitialGapflag = 0; % the first few inputs don't count. Gap it.
    
    global VehicleParams; % for SUV
    VehicleParams.Lf  = 1.12;  % 1.05
    VehicleParams.Lr  = 1.48;  % 1.55
    VehicleParams.L   = 2.6;  %VehiclePara.Lf + VehiclePara.Lr;

    global DataTupleFlag; 
    DataTupleFlag = 0; % the first few inputs don't count. Gap it.  
%  End of mdlInitializeSizes

function sys = mdlUpdates(t,x,u)
%==============================================================
% Update the discrete states, flag = 2， mdlUpdate
% Handle discrete state updates, sample time hits, and major time step
% requirements.
%==============================================================
%  基本没有用到这个过程；在后期的程序模块化时可以继续开发这个功能。
    sys = x;    
% End of mdlUpdate.

function sys = mdlOutputs(t,x,u)
%==============================================================
% Calculate outputs, flag = 3， mdlOutputs
% Return the block outputs. 
%==============================================================

%***********Step (1). Parameters Initialization ***************************************%
% global Initial_State; 
% global P; 
% global PreviousYawrate; 
    
global InitialGapflag;
global VehicleParams;
global params; 
global pre_regressor;
global DataTupleFlag;

Cf_Hat           = 0;
Cr_Hat           = 0;
C_alpha_f_Alt    = 0;
C_alpha_r_Alt    = 0;


t_Elapsed       = 0;
Ay_G_SM         = 0;
Ax              = 0;  
fwa             = 0;
Beta            = 0;

alphaf = 0;
alphar = 0;
Arfa_f = 0;
Arfa_r = 0;

if InitialGapflag < 2 %  get rid of the first two inputs,  because no data from CarSim
    InitialGapflag = InitialGapflag + 1;
else % start control
    InitialGapflag = InitialGapflag + 1;
%***********Step (2). State estimation and Location **********************% 
    t_Start = tic; % 开始计时  
    %-----Update State Estimation of measured Vehicle Configuration--------%
    [Carsim_export] = func_CarsimData_Parse(u);   
    Vx            = Carsim_export.x_dot; 
    Vy            = Carsim_export.y_dot; 
    yawrate       = Carsim_export.phi_dot; % rad/s
    fwa           = Carsim_export.fwa;
    Fyf           = Carsim_export.Fyf;
    Fyr           = Carsim_export.Fyr;   
    alphaf_Direct = Carsim_export.alphaf;
    alphar_Direct = Carsim_export.alphar;
    %-----I. Estimate Cornering stiffness using old method----------------%  
    alphaf_Hat = (Vy + yawrate*VehicleParams.Lf)/Vx - fwa;
    [yf, Calpha_f] = func_RLSFilter_Calpha_f_useFyf(alphaf_Hat, Fyf);
    Cf_Hat = sum(Calpha_f);
    
    %for rear tire 
    alphar_Hat = (Vy - yawrate*VehicleParams.Lr)/Vx;
    [yr, Calpha_r] = func_RLSFilter_Calpha_r_useFyr(alphar_Hat, Fyr);
    Cr_Hat =sum(Calpha_r);

    
    %-----II. Estimate Cornering stiffness using Matrix method----------%  
    if (DataTupleFlag < 10 )
        for J = 10:-1:2
            pre_regressor.input(J)  = pre_regressor.input(J-1);
            pre_regressor.d(J)      = pre_regressor.d(J-1);
        end;
        pre_regressor.input(1)      = alphar_Direct;
        pre_regressor.d(1)          = Fyr;
        DataTupleFlag               = DataTupleFlag + 1;
    else
        for J = 10:-1:2
            pre_regressor.input(J)  = pre_regressor.input(J-1);
            pre_regressor.d(J)      = pre_regressor.d(J-1);
        end;
        pre_regressor.input(1) = alphar_Direct;
        pre_regressor.d(1) = Fyr;
        [outputVector, errorVector, coefficientVector] = func_RLS_Alt_New(pre_regressor.d, pre_regressor.input, params);   
        params.initialCoefficients = coefficientVector;

        C_alpha_f_Alt = 0;
        C_alpha_r_Alt = coefficientVector;    
    end

    
    
    
    
    
    
    t_Elapsed  = toc(t_Start);
     
end % end of if Initialflag < 2 % 

sys = [Cf_Hat; Cr_Hat; C_alpha_f_Alt; C_alpha_r_Alt]; % 

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
% we should do state estimation, but for simplicity we deem that the
% measurements are accurate
% Update the state vector according to the input of the S function,
%           usually do State Estimation from measured Vehicle Configuration
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


