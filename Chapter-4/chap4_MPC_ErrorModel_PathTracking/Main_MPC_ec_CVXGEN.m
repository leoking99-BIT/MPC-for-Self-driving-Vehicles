function [sys,x0,str,ts] =Main_MPC_ec_CVXGEN(t,x,u,flag)
%***************************************************************%
% This is a Simulink/Carsim joint simulation solution for path tracking use
% MPC with tracking error model.
% Use constant high speed, curve path tracking 
% state vector =[epsi,ed,measured_delta_f]
% control input = [steer_SW]
%
% Input:
% t是采样时间, x是状态变量, u是输入(是做成simulink模块的输入,即CarSim的输出),
% flag是仿真过程中的状态标志(以它来判断当前是初始化还是运行等)
%
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
sizes.NumDiscStates  = 2;  %模块离散状态变量的个数,实际上没有用到这个数值，只是用这个来表示离散模块
sizes.NumOutputs     = 11;  %S函数的输出，包括控制量和其它监测量
sizes.NumInputs      = 38; %S函数模块输入变量的个数，即CarSim的输出量
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

%-----------------------------------------------------------------------%
    global InitialGapflag; 
    InitialGapflag = 0; % the first few inputs don't count. Gap it.
    
    global VehiclePara; % for SUV
    VehiclePara.Lf  = 1.05;
    VehiclePara.Lr  = 1.55;
    VehiclePara.L   = 2.6;  %VehiclePara.Lf + VehiclePara.Lr;
%     VehiclePara.Tr  = 1.565;  %c,or 1.57. 注意半轴长度lc还未确定
%     VehiclePara.mu  = 0.85; % 0.55; %地面摩擦因数
%     VehiclePara.Iz  = 2059.2;   %I为车辆绕Z轴的转动惯量，车辆固有参数  
%     VehiclePara.Ix  = 700.7;   %I为车辆绕Z轴的转动惯量，车辆固有参数  
%     VehiclePara.Radius = 0.379;  % 轮胎滚动半径   
    
    global MPCParameters; 
    MPCParameters.Np  = 40;% predictive horizon Assume Np=Nc
    MPCParameters.Ts  = 0.1; % 0.1;  
    MPCParameters.Nx  = 2; %the number of state variables
    MPCParameters.Ny  = 2; %the number of output variables      
    MPCParameters.Nu  = 1; %the number of control inputs
    
    global CostWeights; 
    CostWeights.Wephi   = 5; %state vector =[beta,yawrate,e_phi,s,e_y]
    CostWeights.Wey     = 100;
    CostWeights.deltaf  = 1000;% on Du
    
    global Constraints;  
    Constraints.dumax   = 0.08; %*MPCParameters.Ts; % Units: rad,0.174rad/s = 10deg/s, 0.08rad/s=4.6deg/s  
    Constraints.umax    = 0.471; % Units: rad.  0.4rad=23deg, 0.471rad=27deg
    
    Constraints.DPhimax = pi/3;  %  最大航向角偏差60deg
    Constraints.Dymax   = 1.7; % unit:m. cross-track-error max 2m

    global WayPoints_IndexPre;
    WayPoints_IndexPre = 1;
    
    global Reftraj;
     Reftraj = load('WayPoints_Alt3fromFHWA_Samples.mat');    
    
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


t_Elapsed       = 0;
Ctrl_SteerSW    = 0;
PosX            = 0;
PosY            = 0;
PosPsi          = 0;
Vel             = 0;
e_psi           = 0;
e_d             = 0;
fwa_opt         = 0;
fwa_measured    = 0;
Station         = 0;
    
if InitialGapflag < 2 %  get rid of the first two inputs,  because no data from CarSim
    InitialGapflag = InitialGapflag + 1;
else % start control
    InitialGapflag = InitialGapflag + 1;
%***********Step (2). State estimation and Location **********************% 
    t_Start = tic; % 开始计时  
    %-----Update State Estimation of measured Vehicle Configuration--------%
    [VehStateMeasured, ParaHAT] = func_StateEstimation(u);   
    PosX        = VehStateMeasured.X;
    PosY        = VehStateMeasured.Y;
    PosPsi      = VehStateMeasured.phi;    
    Vel         = VehStateMeasured.x_dot; 
    fwa_measured  = VehStateMeasured.fwa; % rad
    Station     = VehStateMeasured.Station;
    if(Vel < 1.0)
        Vel = 1.0;
    end
    %********Step(3): Given reference trajectory, update vehicle state and bounds *******************% 
    [WPIndex, RefP, RefU, Uaug, PrjP] = func_RefTraj_LocalPlanning( MPCParameters,... 
                            VehiclePara,... 
                            WayPoints_IndexPre,... 
                            Reftraj.WayPoints_Collect,... 
                            VehStateMeasured ); % 
                            
    if ( WPIndex <= 0)
       fprintf('Error: WPIndex <= 0 \n');% 出错
    else
        epsi = PrjP.epsi;  
        if(epsi > pi/2)
           epsi = epsi - pi;
        end
        if(epsi < -pi/2)
           epsi = epsi + pi;
        end
        if(epsi > Constraints.DPhimax)
           epsi = Constraints.DPhimax;
        end
        if(epsi < -Constraints.DPhimax)
           epsi = -Constraints.DPhimax;
        end
        ed   = PrjP.ey;       
        if(ed > Constraints.Dymax)
           ed = Constraints.Dymax;
        end
        if(ed < -Constraints.Dymax)
           ed = -Constraints.Dymax;
        end        
        Xm = [epsi; ed];      
        WayPoints_IndexPre = WPIndex;        
    end

    %****Step(4):  update MPC_error_model_augmented SSM ******************%
    % x(k+1) = Au*x(k)+Bu1*u1 + Bu2 * u2
    [StateSpaceModel] = func_Update_ecMPC_SSM(VehiclePara, MPCParameters, Vel);
    
     %****Step(4):  update Constraints and bounds ********************%
%     Np          = MPCParameters.Np;
%     Eymax       = zeros(Np,1);
%     Eymin       = zeros(Np,1);     
%     LM_right    = -5;
%     LM_middle   = 0;
%     Yroad_L     = -2.5;
%     for i =1:1:Np  % 注意ey是带符号的, Np = 25
%         Eymax(i,1) = (LM_middle - Yroad_L);
%         Eymin(i,1) = (LM_right - Yroad_L);             
%     end
%     [Envelope] = func_ConstraintsBounds(VehiclePara, MPCParameters, Constraints, StateSpaceModel, Vel, CarHat,  Eymax, Eymin); 
    
    %**** Update Cost Weighting Regulation functions ********************%
%     Q = diag([CostWeights.Wephi, CostWeights.Wey]);
%     R = CostWeights.deltaf;
    
    [Q, R] = func_CostWeightingRegulation(CostWeights, Constraints);

    %================CVXGEN solver==================================%
    settings.verbose    = 0;       % 0-Silence; 1-display
    settings.max_iters  = 25;    %Limits the total iterations
    
    params.xm       = Xm;
    params.um       = fwa_measured; % measured front whee angle
    params.uaug     = Uaug;    
    params.Au       = StateSpaceModel.Au;
    params.Bu       = StateSpaceModel.Bu1;
    params.Ba       = StateSpaceModel.Bu2;

    params.Q        = Q;  
    params.R	    = R;     
    
    params.dumax    = Constraints.dumax;
    params.umax     = Constraints.umax;    

    if (40 == MPCParameters.Np)
       [vars, status] = csolve_Np40(params, settings);       
    else
       [vars, status] = csolve_Np20(params, settings);        
    end

    if (1 == status.converged) %if optimization succeeded.
        fwa_opt = vars.u_0; 
%         for i=1:1:20
%             S_opt(i)    = vars.x{i}; 
%             U_opt(i)    = vars.u{i}; 
%         end  
    else
        fwa_opt = vars.u_0;
        fprintf('CVXGEN not-converged at iter= %d\n', InitialGapflag);                  
    end
    
    %====================================================================%
    Ctrl_SteerSW = 19 * fwa_opt*180/pi; % in deg.    
      
    t_Elapsed = toc( t_Start ); %computation time
    
    e_psi           = PrjP.epsi;
    e_d             = PrjP.ey;     

end % end of if Initialflag < 2 % 


sys = [Ctrl_SteerSW; t_Elapsed; PosX; PosY; PosPsi; Station; Vel; e_psi; e_d; fwa_opt; fwa_measured]; %

% end  %End of mdlOutputs.

%==============================================================
% sub functions
%==============================================================    

%***************************************************************%
% **** State estimation
%***************************************************************%
function [VehStatemeasured, HATParameter] = func_StateEstimation(ModelInput)
%***************************************************************%
% we should do state estimation, but for simplicity we deem that the
% measurements are accurate
% Update the state vector according to the input of the S function,
%           usually do State Estimation from measured Vehicle Configuration
%***************************************************************%  
    %******输入接口转换***%        
    g = 9.81;
    VehStatemeasured.X       = round(100*ModelInput(1))/100;%单位为m, 保留2位小数
    VehStatemeasured.Y       = round(100*ModelInput(2))/100;%单位为m, 保留2位小数    
    VehStatemeasured.phi     = (round(10*ModelInput(3))/10)*pi/180; %航向角，Unit：deg-->rad，保留1位小数    
    VehStatemeasured.x_dot   = ModelInput(4)/3.6; %Unit:km/h-->m/s，保留1位小数  
    VehStatemeasured.y_dot   = ModelInput(5)/3.6; %Unit:km/h-->m/s，保留1位小数   
    VehStatemeasured.phi_dot = (round(10*ModelInput(6))/10)*pi/180; %Unit：deg/s-->rad/s，保留1位小数      
    VehStatemeasured.beta    = (round(10*ModelInput(7))/10)*pi/180;% side slip, Unit:deg-->rad，保留1位小数    
    VehStatemeasured.delta_f = (round(10*0.5*(ModelInput(8)+ ModelInput(9)))/10); % deg
    VehStatemeasured.fwa     = VehStatemeasured.delta_f * pi/180;  % deg-->rad
    VehStatemeasured.Steer_SW= ModelInput(10); %deg
    VehStatemeasured.Ax      = g*ModelInput(11);%单位为m/s^2, 保留2位小数
    VehStatemeasured.Ay      = g*ModelInput(12);%单位为m/s^2, 保留2位小数
    VehStatemeasured.yawrate_dot = ModelInput(13); %rad/s^2
    % Here I don't explore the state estimation process, and deem the
    % measured values are accurate!!! 
    HATParameter.alpha_l1   = (round(10*ModelInput(14))/10)*pi/180; % deg-->rad，保留1位小数   
    HATParameter.alpha_l2   = (round(10*ModelInput(15))/10)*pi/180; % deg-->rad，保留1位小数   
    HATParameter.alpha_r1   = (round(10*ModelInput(16))/10)*pi/180; % deg-->rad，保留1位小数   
    HATParameter.alpha_r2   = (round(10*ModelInput(17))/10)*pi/180; % deg-->rad，保留1位小数     
    HATParameter.alphaf     = (round(10*0.5 * (ModelInput(14)+ ModelInput(16)))/10)*pi/180; % deg-->rad，保留1位小数   
    HATParameter.alphar     = (round(10*0.5 * (ModelInput(15)+ ModelInput(17)))/10)*pi/180; % deg-->rad，保留1位小数  
    
    HATParameter.Fz_l1      = round(10*ModelInput(18))/10; % N 
    HATParameter.Fz_l2      = round(10*ModelInput(19))/10; % N 
    HATParameter.Fz_r1      = round(10*ModelInput(20))/10; % N 
    HATParameter.Fz_r2      = round(10*ModelInput(21))/10; % N 
    
    HATParameter.Fy_l1      = round(10*ModelInput(22))/10; % N 
    HATParameter.Fy_l2      = round(10*ModelInput(23))/10; % N 
    HATParameter.Fy_r1      = round(10*ModelInput(24))/10; % N 
    HATParameter.Fy_r2      = round(10*ModelInput(25))/10; % N 
    HATParameter.Fyf        = HATParameter.Fy_l1 + HATParameter.Fy_r1;
    HATParameter.Fyr        = HATParameter.Fy_l2 + HATParameter.Fy_r2;
    
    HATParameter.Fx_L1      = ModelInput(26);
    HATParameter.Fx_L2      = ModelInput(27);
    HATParameter.Fx_R1      = ModelInput(28);
    HATParameter.Fx_R2      = ModelInput(29);
    
%     HATParameter.GearStat    = ModelInput(30);
    VehStatemeasured.Roll_Shad   = ModelInput(30)*pi/180;% deg-->rad 
    HATParameter.Roll        = ModelInput(31)*pi/180;% deg-->rad 
    HATParameter.Rollrate    = ModelInput(32)*pi/180;% deg/s-->rad/s
    HATParameter.Roll_accel  = ModelInput(33); % rad/s^2
    HATParameter.Z0          = ModelInput(34); %m
    VehStatemeasured.Station     = ModelInput(35); %m
    HATParameter.Zcg_TM      = ModelInput(35); %m
    HATParameter.Zcg_SM      = ModelInput(36); %m
    HATParameter.Ay_CG       = ModelInput(37)*g; %m/s^2
    HATParameter.Ay_Bf_SM    = ModelInput(38)*g; %m/s^2
    
% end % end of func_StateEstimation


