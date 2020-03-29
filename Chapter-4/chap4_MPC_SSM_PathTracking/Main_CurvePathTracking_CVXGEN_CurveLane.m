function [sys,x0,str,ts] =Main_CurvePathTracking_CVXGEN_CurveLane(t,x,u,flag)
%***************************************************************%
% This is a Simulink/Carsim joint simulation solution for safe driving
% envelope control of high speed autonomous vehicle
% Linearized spatial bicycle vehicle dynamic model is applied.
% No successive linearizarion. No two time scale of prediction horizon
% Constant high speed, curve path tracking 
% state vector =[beta,yawrate,e_phi,s,e_y]
% control input = [steer_SW]
% many other parameters are also outputed for comparision.

% Input:
% t是采样时间, x是状态变量, u是输入(是做成simulink模块的输入,即CarSim的输出),
% flag是仿真过程中的状态标志(以它来判断当前是初始化还是运行等)

% Output:
% sys输出根据flag的不同而不同(下面将结合flag来讲sys的含义), 
% x0是状态变量的初始值, 
% str是保留参数,设置为空
% ts是一个1×2的向量, ts(1)是采样周期, ts(2)是偏移量
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
% Return the sizes, initial conditions, and sample times for the S-function.
%==============================================================
sizes = simsizes;%用于设置模块参数的结构体用simsizes来生成
sizes.NumContStates  = 0;  %模块连续状态变量的个数
sizes.NumDiscStates  = 3;  %模块离散状态变量的个数,实际上没有用到这个数值，只是用这个来表示离散模块
sizes.NumOutputs     = 8;  %S函数的输出，包括控制量和其它监测量
sizes.NumInputs      = 7; %S函数模块输入变量的个数，即CarSim的输出量
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

%--Global parameters and initialization
    global InitialGapflag; 
    InitialGapflag = 0; % the first few inputs don't count. Gap it.
    
    global VehiclePara; 
    VehiclePara.m  = 1540;   %m为车辆质量,Kg; Sprung mass = 1370
    VehiclePara.g  = 9.8;
    VehiclePara.Lf = 1.11;  %a
    VehiclePara.Lr = 1.67;  %b，前后车轮距离车辆质心的距离，车辆固有参数
    VehiclePara.L  = 2.78;  %VehiclePara.Lf + VehiclePara.Lr;
    VehiclePara.Lc = 1.59;  %c,or 1.57. 注意半轴长度lc还未确定
    VehiclePara.I  = 2315.3;   %I为车辆绕Z轴的转动惯量，车辆固有参数
    VehiclePara.mu = 1.0; % 0.55; %地面摩擦因数，
    VehiclePara.Radius = 0.261;  % 轮胎滚动半径
    
    global MPCParameters; 
    MPCParameters.Np  = 40;% predictive horizon Assume Np=Nc
    MPCParameters.Ts  = 0.05; %Set the sample time of near term 
    MPCParameters.Nx  = 3; %the number of state variables
    MPCParameters.Ny  = 3; %the number of output variables      
    MPCParameters.Nu  = 2; %the number of control inputs
    
    global CostWeights; 
    CostWeights.Q   = [ 10      0       0;
                        0      10       0;
                        0      0       10];  %state vector =[beta,yawrate,e_phi,s,e_y]
 
    CostWeights.R   = 10000; % on Du
    
    global Constraints;  
    Constraints.dumax   = 0.08;      % Units: rad  
    Constraints.umax    = 0.4;      % Units: rad appro.23deg
    
    global WayPoints_IndexPre;
    WayPoints_IndexPre = 1;
    
    global Reftraj;
    Reftraj = load('WayPoints_Alt3fromFHWA_Overall.mat');
%     Local_reftraj = load('WayPoints_Alt3fromFHWA_Portion.mat');    
    
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

global InitialGapflag;
global VehiclePara;
global MPCParameters; 
global CostWeights;     
global Constraints;
global WayPoints_IndexPre;
global Reftraj;

Ctrl_SteerSW    = 0;
t_Elapsed       = 0;
PosX            = 0;
PosY            = 0;
PosPsi          = 0;
Vel             = 0;
e_psi           = 0;
e_y             = 0;

if InitialGapflag < 2 %  get rid of the first two inputs,  because no data from CarSim
    InitialGapflag = InitialGapflag + 1;
else % start control
%***********Step (2). State estimation and Location **********************% 
    %-----Update State Estimation of measured Vehicle Configuration--------%
    [VehStateMeasured] = func_StateEstimation(u); %u是S函数模块的输入
    Vel     = VehStateMeasured.Vx;     
    PosX    = VehStateMeasured.Xc;
    PosY    = VehStateMeasured.Yc;
    PosPsi  = VehStateMeasured.psi;    
    Ax      = VehStateMeasured.Ax;
    fwa     = VehStateMeasured.fwa;    

    %********Step(3): Given reference trajectory, update vehicle state and bounds *******************% 
    % Local vehicle coordinates, y=K[3]x^3 + K[2]x^2 + K[1]x +K[0].
    [PrjP, RefP, RefU, WPIndex] = func_RefTraj_LocalPlanning( MPCParameters, WayPoints_IndexPre, Reftraj.WayPoints_Collect, VehStateMeasured ); % reference path is a straight line
    if ( WPIndex <= 0)
        %出错
    else
        WayPoints_IndexPre = WPIndex;        
    end

    %****Step(4):  MPC formulation;********************%
    %----Update  An, Al, B, dn,dl of the StateSpaceModel
    [StateSpaceModel] = func_StateSpaceModel_StraightLane(VehiclePara, MPCParameters,  PrjP ); 
    Xm = [PosX; PosY; PosPsi];
    
    %================CVXGEN solver==================================%
    settings.verbose    = 0;       % 0-Silence; 1-display
    settings.max_iters  = 25;    %Limits the total iterations
    
    params.xm       = Xm;
    params.um       = fwa; % measured front whee angle
    params.Pxr      = [PrjP.xr; PrjP.yr; PrjP.psir];
    params.Pur      = PrjP.fwar;
    params.An       = StateSpaceModel.An;
    params.Bn       = StateSpaceModel.Bn;
    params.Q        = CostWeights.Q;  
    params.R	    = CostWeights.R;
    params.umax     = Constraints.umax;
    params.dumax    = Constraints.dumax; 
    params.RefP     = RefP; 
    params.RefU     = RefU;  
   
    t_Start = tic; % 开始计时  
    [vars, status] = csolve_StraightLane(params, settings);
    if (1 == status.converged) %if optimization succeeded.
        fwa_opt = vars.u_0;          
%         ah_des  = vars.u_0(2); 
    else
        fwa_opt =  0;
        fprintf('CVXGEN converged = 0 \n');                  
    end
    
    %====================================================================%
    Ctrl_SteerSW0 = 19 * fwa_opt*180/pi; % in deg.    
%     [Throttle, Brake] = func_AccelerationTrackingController(ah_opt);
      
    t_Elapsed = toc( t_Start ); %computation time
    
     %---4.Publish command********************%
    Ctrl_SteerSW = round(10*Ctrl_SteerSW0)/10;
    e_y            = PrjP.ey;
    e_psi          = PrjP.epsi;   
end % end of if Initialflag < 2 % 

   
sys = [Ctrl_SteerSW; t_Elapsed; PosX; PosY; PosPsi; Vel; e_psi; e_y];     
% end  %End of mdlOutputs.

%==============================================================
% sub functions
%==============================================================    

%***************************************************************%
% **** State estimation
%***************************************************************%
function [VehStatemeasured] = func_StateEstimation(ModelInput)
%***************************************************************%
% we should do state estimation, but for simplicity we deem that the
% measurements are accurate
% Update the state vector according to the input of the S function,
%           usually do State Estimation from measured Vehicle Configuration
%***************************************************************%  
    %******输入接口转换***%        
    VehStatemeasured.Vx      = ModelInput(1)/3.6; %Unit:km/h-->m/s，保留1位小数  
    VehStatemeasured.Xc      = round(100*ModelInput(2))/100;%单位为m, 保留2位小数
    VehStatemeasured.Yc      = round(100*ModelInput(3))/100;%单位为m, 保留2位小数    
    VehStatemeasured.psi     = (round(100*ModelInput(4))/100)*pi/180; %航向角，Unit：deg-->rad，保留2位小数    
    VehStatemeasured.Ax      = 9.8*ModelInput(5);%单位为m/s^2, 保留2位小数
    VehStatemeasured.fwa     = (round(10*0.5*(ModelInput(6)+ ModelInput(7)))/10)*pi/180; % deg-->rad   
    
% end % end of func_StateEstimation

%***************************************************************%
% Augmented vehicle state space model
%***************************************************************%
function [StateSpaceModel] = func_StateSpaceModel_StraightLane(VehiclePara, MPCParameters, PrjP)
    % generate State-space model
    L       = VehiclePara.L;  %a = 1.11;  
    Ts      = MPCParameters.Ts;
    Nx      = MPCParameters.Nx;
    Velr    = PrjP.Velr;    
    xr      = PrjP.xr;
    yr      = PrjP.yr;
    psir    = PrjP.psir;  
    fwar    = PrjP.fwar; 
   
    
    Acn = [0,        0,        -Velr*sin(psir);
           0,        0,        Velr*cos(psir);
           0,        0,        0];
      
    Bcn = [0,       0,     Velr/(L*cos(fwar)*cos(fwar))]';%
   
    % SSM discretization for the near term
    Adn = eye(Nx) + Ts*Acn;
    Bdn = Ts*Bcn;
 
    StateSpaceModel.An = Adn;
    StateSpaceModel.Bn = Bdn; 
                        
% end % end of func_SpatialDynamicalModel


