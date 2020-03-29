function [sys,x0,str,ts] =Main_KinematicModel_Validation(t,x,u,flag)
%***************************************************************%
% 第2.1节所建立的车辆运动模型的仿真验证
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
% end %  end sfuntmpl

%==============================================================
% Initialization, flag = 0，mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function.
%==============================================================
function [sys,x0,str,ts] = mdlInitializeSizes
%***************************************************************%
% Call simsizes for a sizes structure, fill it in, and convert it 
% to a sizes array.
%***************************************************************% 
sizes = simsizes;%用于设置模块参数的结构体用simsizes来生成
sizes.NumContStates  = 0;  %模块连续状态变量的个数
sizes.NumDiscStates  = 5;  %模块离散状态变量的个数,实际上本app 没有用到这个数值，只是用这个来表示离散模块
sizes.NumOutputs     = 10;  %S函数的输出，包括控制量和其它监测量
sizes.NumInputs      = 10; %S函数模块输入变量的个数，即CarSim的输出量
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
  
global InitialGapflag; 
    InitialGapflag = 0; % Ignore the first few inputs from CarSim
global Previous_States; % store the previous state vector 
    Previous_States.X_pred = 0.0; 
    Previous_States.Y_pred = 0.0; 
    Previous_States.Yaw_pred = 0.0;

% RLS initialization
[y, c] = func_SteerRatio_Estimation_RLS('initial', 0.95, 1, 10);
[y, e] = func_SteerRatio_Estimation_RLS_array('initial', 0.95, 10, 10);
%  End of mdlInitializeSizes

%==============================================================
% Update the discrete states, flag = 2， mdlUpdate
% Handle discrete state updates, sample time hits, and major time step
% requirements.
%==============================================================
function sys = mdlUpdates(t,x,u)
%  目前没有用到这个过程；在后期的程序模块化时可以继续开发这个功能。
    sys = x;
% end     %End of mdlUpdate.

%==============================================================
% Calculate outputs, flag = 3， mdlOutputs
% Return the block outputs. 
%==============================================================
function sys = mdlOutputs(t,x,u)
%t是采样时间, x是状态变量, u是输入(是做成simulink模块的输入)
global InitialGapflag; 
global Previous_States;
lfr = 2.78;
Ts = 0.05;
Steer_ratio = 1;

    % 提取CarSim输入到Simulink的数据
    x_L2 = u(1); %左后轮x坐标
    x_R2 = u(2); %右后轮x坐标
    y_L2 = u(3); %左后轮y坐标
    y_R2 = u(4); %右后轮y坐标   
    Yaw  = u(5)*pi/180;%航向角Unit：deg-->rad
    Steer_SW = u(6); %方向盘角度
    Steer_L1 = u(7); %左前轮偏角
    Steer_R1 = u(8); %右前轮偏角
    Vx_L2 = u(9);  %左后轮纵向速度，Unit:km/h
    Vx_R2 = u(10); %右后轮纵向速度，Unit:km/h
    
    Car_X = 0.5*(x_L2 + x_R2);%后轴中心X坐标，Unit:m
    Car_Y = 0.5*(y_L2 + y_R2);%后轴中心Y坐标，Unit:m
    Vx_km_h = 0.5*(Vx_L2 + Vx_R2);%后轴中心处纵向速度,Unit：km/h
    Steer_deg = 0.5*(Steer_L1 + Steer_R1);%等效前轮偏角，Unit：deg
    
    Vx_m_s  = Vx_km_h/3.6;%%后轴中心处纵向速度 in (m/s),Unit：m/s    
    Steer_rad = Steer_deg*pi/180;%等效前轮偏角in (rad)，Unit：degs-->rad;

if (InitialGapflag < 3) %  Ignore the first few inputs
    InitialGapflag = InitialGapflag + 1;
    X_pred = Car_X; 
    Y_pred = Car_Y; 
    Yaw_pred = Yaw;
    Previous_States.X_pred = Car_X; 
    Previous_States.Y_pred = Car_Y; 
    Previous_States.Yaw_pred = Yaw;
else % start control
    %-----I. Update predicted states using differential equation--------%
%     Updated_state = func_UpdateState_dsolve_2_7(Previous_States, lfr, Vx_m_s, Steer_rad, Ts);

    %-----II. Update predicted states using RK4-------%
    
    %-----III. Update predicted states using Euler Method------%
    Updated_state = func_UpdateState_EulerM_2_7(Previous_States, lfr, Vx_m_s, Steer_rad, Ts);
    
    X_pred = Updated_state.X_pred; 
    Y_pred = Updated_state.Y_pred; 
    Yaw_pred = Updated_state.Yaw_pred;
    
    Previous_States.X_pred = X_pred; 
    Previous_States.Y_pred = Y_pred; 
    Previous_States.Yaw_pred = Yaw_pred;

    %-----Estimate Steer_ratio-----%  
%     [Steer_SW_hat, Steer_ratio, Rinv_f] = func_SteerRatio_Estimation_RLS(Steer_deg, Steer_SW);
%     Hat_err = Steer_SW_hat - Steer_SW;
    [Steer_SW_hat, Steer_ratio_vector] = func_SteerRatio_Estimation_RLS_array(Steer_deg, Steer_SW);
    Steer_ratio =sum(Steer_ratio_vector);
 
 
end % End of if (Initialflag < 3) % 

    
    sys = [Car_X; Car_Y; Yaw; X_pred; Y_pred; Yaw_pred; Vx_m_s; Steer_rad; Steer_SW; Steer_ratio];  
       
% end  %End of mdlOutputs.

