function [sys,x0,str,ts] =Main_MPC_SpeedCtrl(t,x,u,flag)
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

%==============================================================
% Initialization, flag = 0，mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function.
%==============================================================
function [sys,x0,str,ts] = mdlInitializeSizes
sizes = simsizes;%用于设置模块参数的结构体用simsizes来生成
sizes.NumContStates  = 0;  %模块连续状态变量的个数
sizes.NumDiscStates  = 2;  %模块离散状态变量的个数,实际上没有用到这个数值，只是用这个来表示离散模块
sizes.NumOutputs     = 6;  %S函数的输出，包括控制量和其它监测量
sizes.NumInputs      = 2; %S函数模块输入变量的个数，即CarSim的输出量
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
global MPCParameters; 
    MPCParameters.Np      = 30;% predictive horizon
    MPCParameters.Nc      = 30;% control horizon
    MPCParameters.Nx      = 2; %number of state variables
    MPCParameters.Nu      = 1; %number of control inputs
    MPCParameters.Ny      = 1; %number of output variables  
    MPCParameters.Ts      = 0.05; %Set the sample time
    MPCParameters.Q       = 100; % cost weight factor 
    MPCParameters.R       = 0.1; % cost weight factor 
    MPCParameters.S       = 0.1; % cost weight factor 
    MPCParameters.qp_solver = 0; %0: default, quadprog; 1:qpOASES
    MPCParameters.refspeedT = 1; %0: default, step speed profile; 
                                 %1:sine-wave speed profile
    MPCParameters.umin      = -5.0;  % the min of deceleration
    MPCParameters.umax      = 3.5;  % the max of acceleration
    MPCParameters.dumin     = -5.0; % minimum limits of jerk
    MPCParameters.dumax     = 5.0; % maximum limits of jerk
global WarmStart;
    WarmStart = zeros(MPCParameters.Np,1);
%  End of mdlInitializeSizes

%==============================================================
% Update the discrete states, flag = 2， mdlUpdate
% Handle discrete state updates, sample time hits, and major time step
% requirements.
%==============================================================
function sys = mdlUpdates(t,x,u)
%  基本没有用到这个过程；在后期的程序模块化时可以继续开发这个功能。
    sys = x;    
% End of mdlUpdate.

%==============================================================
% Calculate outputs, flag = 3， mdlOutputs
% Return the block outputs. 
%==============================================================
function sys = mdlOutputs(t,x,u)

global InitialGapflag;
global MPCParameters;
global WarmStart;
Vx    = 0;
a_x   = 0;
a_des = 0;

t_Start = tic; % 开始计时 
if InitialGapflag < 2 %  get rid of the first two inputs
    InitialGapflag = InitialGapflag + 1;% delay
else
    InitialGapflag = InitialGapflag + 1;
    %***********Step (1). Update vehicle states *************************% 
    Vx    = u(1)/3.6;  %车辆纵向速度，单位：km/h-->m/s
    a_x   = u(2)*9.8;  %车辆纵向加速度，单位：g's-->m/s2 
    kesi = [Vx;  a_x]; %更新车辆状态向量
    
    %********Step(2): Generate reference speed profile *******************%
    switch MPCParameters.refspeedT,
        case 0 % default, step speed profile
            %----设定阶梯式的期望速度曲线----------------------%
            SpeedProfile = func_ConstantSpeed(InitialGapflag, MPCParameters);
        case 1 % sine-wave speed profile
            %----设置sine wave形式 的期望速度曲线--------------%
            SpeedProfile = func_SineSpeed(InitialGapflag,MPCParameters);
        otherwise % Unexpected flags %
            error(['unexpected speed-profile:',num2str(MPCParameters.refspeedT)]); % Error handling
    end %  end of switch 
    
    %****Step(3): update longitudinal vehilce model with inertial delay**8%
    Ts = MPCParameters.Ts; % 50ms
    tao = 0.2; %车辆纵向系统惯性延时参数，0.2 for simulation
    StateSpaceModel.A =  [1      Ts;
                          0      1-Ts/tao];
    StateSpaceModel.B =  [0;     Ts/tao]; 
    StateSpaceModel.C =  [1,    0];
    %****Step(4):  MPC formulation;********************%
    %Update Theta and PHI for future states prediction
    [PHI, THETA] = func_Update_PHI_THETA(StateSpaceModel, MPCParameters);

    %Update H and f for cost function J
    [H, f, g] = func_Update_H_f(kesi, SpeedProfile, PHI, THETA, MPCParameters); 
    
    %****Step(5):  Call qp-solver********************%
    switch MPCParameters.qp_solver,
        case 0 % default qp-solver: quadprog
            [A, b, Aeq, beq, lb, ub] = func_Constraints_u_quadprog(MPCParameters);
            options = optimset('Display','off', ...
                            'TolFun', 1e-8, ...
                            'MaxIter', 2000, ...
                            'Algorithm', 'active-set', ...
                            'FinDiffType', 'forward', ...
                            'RelLineSrchBnd', [], ...
                            'RelLineSrchBndDuration', 1, ...
                            'TolConSQP', 1e-8); 
            warning off all  % close the warnings during computation     

            U0 = WarmStart;           
            [U, FVAL, EXITFLAG] = quadprog(H, g, A, b, Aeq, beq, lb, ub, U0, options); %
            WarmStart = shiftHorizon(U);     % Prepare restart, nominal close loop 
            if (1 ~= EXITFLAG) %if optimization NOT succeeded.
                U(1) = 0.0;
                fprintf('MPC solver not converged!\n');                  
            end
            a_des =  U(1);
 
        case 1 % qpOASES
            [A, lb, ub, lbA, ubA] = func_Constraints_u_qpOASES(MPCParameters);
            options = qpOASES_options('default', ...
                                'printLevel', 0); 

            %=======================USE QP==================%
            [U, FVAL, EXITFLAG, iter, lambda] = qpOASES(H, g, A, lb, ub, lbA, ubA, options); %

            %=======================USE SQP==================%
    %         try
    %             H=sparse(H);
    %             A=sparse(A);
    %         catch
    %             fprintf('qpOASES Error reported\n'); 
    %         end
    %         if (qpOASES_hotstart_flag)
    %             [qpOASES_QP, U, FVAL, EXITFLAG, iter, lambda] = qpOASES_sequence('i', H, g, A, lb, ub, lbA, ubA, options);
    %             qpOASES_hotstart_flag = 1;
    %         else    
    %             [U, FVAL, EXITFLAG, iter, lambda] = qpOASES_sequence('m', qpOASES_QP, H, g, A, lb, ub, lbA, ubA, options); %
    %         end
            if (0 ~= EXITFLAG) %if optimization NOT succeeded.
                U(1) = 0.0;
                fprintf('MPC solver: qpOASES not converged!\n');                  
            end
            a_des =  U(1);

        otherwise % Unexpected flags %
            error(['unexpected qp-solver, Sol_method=',num2str(flag)]); % Error handling
    end %  end of switch 

end % end of if Initialflag < 1 % 

    %****Step(6):  由期望的加速度生成Throttle和Brake;********************%
    [Throttle, Brake] = func_AccelerationTrackingController(a_des);

t_Elapsed = toc( t_Start ); %computation time 

sys = [Throttle; Brake;t_Elapsed; Vx; a_x; a_des]; 
% end  %End of mdlOutputs.

%==============================================================
% sub functions
%============================================================== 
function [Vref] = func_SineSpeed(Index, MPCParameters)
%生成正弦形式的期望速度曲线
    %****Sine wave parameters
    T = 50; %正弦速度曲线的周期，unit: s
    freq = 1/T; %正弦速度曲线的频率，unit: Hz
    Amplit = 10;%正弦速度曲线的幅值
    offst = 20; %正弦速度曲线的偏移
    
    Ts = MPCParameters.Ts; %采样时间=0.05，unit: s
    Np = MPCParameters.Np; % 预测时域：30
    Vref = cell(Np,1);
    t0 = Index*Ts;

    for i = 1:1:Np
        t = t0 + i*Ts;
        Vref{i,1}   =   Amplit*sin(2*pi*freq*t) + offst;   
    end
   
% end %EoF

function [Vref] = func_ConstantSpeed(InitialGapflag, MPCParameters)
% 生成阶梯形式的期望速度曲线    
    Ts = MPCParameters.Ts; %采样时间=0.05，unit: s
    Np = MPCParameters.Np; % 预测时域：30
    Vref = cell(Np,1);
    
    % 自定义阶梯的形式
    if InitialGapflag < 400
        Vset = 10;
    else
        if InitialGapflag < 800
            Vset = 10;
        else
            if InitialGapflag < 1500
                Vset = 20;
            else
                Vset = 5;
            end
        end
    end

    for i = 1:1:Np
        Vref{i,1}   =   Vset;   
    end

% end %EoF

function [Throttle, Brake] = func_AccelerationTrackingController(ahopt)
% 车辆下位控制器将期望加速度转化为油门控制量和制动主缸压力控制量
    K_brake         = 0.3;
    K_throttle      = 0.1; %0.05;
    Brake_Sat       = 15;
    Throttle_Sat    = 1;

    if ahopt < 0 % Brake control
        Brake = K_brake * ahopt;
        if Brake > Brake_Sat
            Brake = Brake_Sat;
        end
        Throttle = 0;
    else % throttle control 
        Brake       = 0;
        Throttle    = K_throttle  *ahopt;
        if Throttle > Throttle_Sat
            Throttle = Throttle_Sat;
        end
        if Throttle < 0
            Throttle = 0;
        end

    end
% end %EoF

function u0 = shiftHorizon(u) %shift control horizon
    u0 = [u(:,2:size(u,2)), u(:,size(u,2))];  %  size(u,2))
    
function [PHI, THETA] = func_Update_PHI_THETA(StateSpaceModel, MPCParameters)
%***************************************************************%
% 预测输出表达式 Y(t)=PHI*kesi(t)+THETA*DU(t) 
% Y(t) = [Eta(t+1|t) Eta(t+2|t) Eta(t+3|t) ... Eta(t+Np|t)]'
%***************************************************************%
    Np = MPCParameters.Np;
    Nc = MPCParameters.Nc;
    Nx = MPCParameters.Nx;
    Ny = MPCParameters.Ny;
    Nu = MPCParameters.Nu;
    A = StateSpaceModel.A;
    B = StateSpaceModel.B;
    C = StateSpaceModel.C;

    PHI_cell=cell(Np,1);                            %PHI=[CA CA^2  CA^3 ... CA^Np]' 
    THETA_cell=cell(Np,Nc);                         %THETA
    for j=1:1:Np
        PHI_cell{j,1}=C*A^j;                       %  demision:Ny* Nx
        for k=1:1:Nc
            if k<=j
                THETA_cell{j,k}=C*A^(j-k)*B;        %  demision:Ny*Nu
            else 
                THETA_cell{j,k}=zeros(Ny,Nu);
            end
        end
    end
    PHI=cell2mat(PHI_cell);    % size(PHI)=[(Ny*Np) * Nx]
    THETA=cell2mat(THETA_cell);% size(THETA)=[Ny*Np Nu*Nc]
% end %EoF


function[H, f, g] = func_Update_H_f(kesi, SpeedProfile, PHI, THETA, MPCParameters)
%***************************************************************%
% trajectory planning
%***************************************************************%
    Np = MPCParameters.Np;
    Nc = MPCParameters.Nc;   
    Q  = MPCParameters.Q;
    R  = MPCParameters.R;
        
    Qq = kron(eye(Np),Q);  %           Q = [Np*Nx] *  [Np*Nx] 
    Rr = kron(eye(Nc),R);  %           R = [Nc*Nu] *  [Nc*Nu]

    Vref = cell2mat(SpeedProfile); 
    error = PHI * kesi;    %[(Nx*Np) * 1]

    H = THETA'*Qq*THETA + Rr;  
    f = (error' - Vref')*Qq*THETA;
    g = f';
% end %EoF

function  [A, b, Aeq, beq, lb, ub] = func_Constraints_u_quadprog(MPCParameters)
%************************************************************************%
% generate the constraints of the vehicle
%  
%************************************************************************%
    Np   = MPCParameters.Np;
    Nc   = Np;    
    umin = MPCParameters.umin;
    umax = MPCParameters.umax;  

    A   = [];
    b   = [];
    Aeq = [];
    beq = [];

%----(3) lb=<x<=ub----------%
    lb=kron(ones(Nc,1),umin);
    ub=kron(ones(Nc,1),umax);
% end %EoF

function [A_t, lb, ub, lbA, ubA] = func_Constraints_u_qpOASES(MPCParameters)
    Np   = MPCParameters.Np;
    Nc   = Np;    
    umin = MPCParameters.umin;
    umax = MPCParameters.umax;  

    A_t = [];
    lbA = [];
    ubA = [];
%---- lb=<x<=ub----------%
    lb=kron(ones(Nc,1),umin);
    ub=kron(ones(Nc,1),umax);
% end %EoF
