function [sys,x0,str,ts] =Main_MPC_Dynamics_AFS(t,x,u,flag)
%***************************************************************%
% 该程序功能：用线性化的车辆动力学模型（小角度假设）设计MPC控制器，
% 实现不同路面情况、车速下的双移线跟踪，并利用Simulink/CarSim实现联合仿真
% MATLAB版本：R2013b,CarSim版本：8.1
% 状态量=[y_dot,x_dot,phi,phi_dot,Y,X]，控制量为前轮偏角delta_f

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
sizes.NumDiscStates  = 6;  %模块离散状态变量的个数,实际上没有用到这个数值，只是用这个来表示离散模块
sizes.NumOutputs     = 7;  %S函数的输出，包括控制量和其它监测量
sizes.NumInputs      = 25; %S函数模块输入变量的个数，即CarSim的输出量
sizes.DirFeedthrough = 1;  %模块是否存在直接贯通(direct feedthrough). 
sizes.NumSampleTimes = 1;  %模块的采样次数，>=1
sys = simsizes(sizes);    %设置完后赋给sys输出
x0 = zeros(sizes.NumDiscStates,1);%initial the  state vector， of no use

str = [];             % 保留参数，Set str to an empty matrix.
ts  = [0.05 0];       % ts=[period, offset].采样周期sample time=0.05s 

%--Global parameters and initialization
    % [y, e] = func_RLSFilter_Ccf('initial', 0.95, 10, 10);
    % [y, e] = func_RLSFilter_Ccr('initial', 0.95, 10, 10);
    % [y, e] = func_RLSFilter_Clf('initial', 0.95, 10, 10);
    % [y, e] = func_RLSFilter_Clr('initial', 0.95, 10, 10);

global InitialGapflag; 
    InitialGapflag = 0; % the first few inputs don't count. Gap it.
    
global VehiclePara; % for SUV
    VehiclePara.m   = 1600;   %m为车辆质量,Kg; Sprung mass = 1370
    VehiclePara.g   = 9.81;
    VehiclePara.hCG = 0.65;%m
    VehiclePara.Lf  = 1.12;  % 1.05
    VehiclePara.Lr  = 1.48;  % 1.55
    VehiclePara.L   = 2.6;  %VehiclePara.Lf + VehiclePara.Lr;
    VehiclePara.Tr  = 1.565;  %c,or 1.57. 注意半轴长度lc还未确定
    VehiclePara.mu  = 0.85; % 0.55; %地面摩擦因数
    VehiclePara.Iz  = 2059.2;   %I为车辆绕Z轴的转动惯量，车辆固有参数  
    VehiclePara.Ix  = 700.7;   %I为车辆绕Z轴的转动惯量，车辆固有参数  
    VehiclePara.Radius = 0.387;  % 轮胎滚动半径   
    
global MPCParameters; 
    MPCParameters.Np  = 25;% predictive horizon Assume Np=Nc
    MPCParameters.Nc  = 15; %  Tsplit
    MPCParameters.Ts  = 0.05; % the sample time of near term  
    MPCParameters.Nx  = 6; %the number of state variables
    MPCParameters.Ny  = 2; %the number of output variables      
    MPCParameters.Nu  = 1; %the number of control inputs
    
global WarmStart;
    WarmStart = zeros(MPCParameters.Nu * MPCParameters.Nc,1);
    
global CostWeights; 
    CostWeights.Wephi    = 100; %state vector =[beta,yawrate,e_phi,s,e_y]
    CostWeights.Wey      = 100;
    CostWeights.WDdeltaf = 1000;

global Constraints;
%     Constraints.dumax   = 0.08; % Units: rad,0.08rad=4.6deg
%     Constraints.dumax   = 0.1; % Units: rad,0.1rad=5.7deg  
    Constraints.dumax   = 0.0148; % Units: rad,0.0148rad = 0.8deg
    Constraints.umax    = 0.4;  % Units: rad, 0.4rad=23deg
    
    Constraints.ycmin   = [-0.5;  -5];
    Constraints.ycmax   = [0.5;   5];

global WayPoints_IndexPre;
    WayPoints_IndexPre = 1;

global Reftraj;
    Reftraj = load('Waypoints_Double_Line_Shift.mat'); 
    
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
global WarmStart;
global CostWeights;     
global Constraints;
global WayPoints_IndexPre;
global Reftraj;
    
Ts = MPCParameters.Ts;
Np = MPCParameters.Np;
Nc = MPCParameters.Nc;
Nx = MPCParameters.Nx;
Ny = MPCParameters.Ny;
Nu = MPCParameters.Nu;
Naug = Nx + Nu;

Steer_SW_deg    = 0;
t_Elapsed       = 0;
PosX            = 0;
PosY            = 0;
PosPhi          = 0;
e_psi           = 0;
e_y             = 0;

if InitialGapflag < 2 %  get rid of the first two inputs
    InitialGapflag = InitialGapflag + 1;
else % start control
    InitialGapflag = InitialGapflag + 1;
    %**********Step (2). Update state and tire-stiffness estimation ******% 
    t_Start = tic; % 开始计时  
    %-----Update State Estimation of measured Vehicle Configuration-------%
    y_dot = u(1)/3.6; %CarSim输出的是km/h，转换为m/s
    x_dot = u(2)/3.6;%CarSim输出的是km/h，转换为m/s
    if (0 == x_dot) %若x_dot=0，则将其设为一个非常小的数，防止出现分母为零的情况
        x_dot = 0.001;
    end
    PosPhi = u(3)*pi/180; %CarSim输出的为角度，角度转换为弧度
    phi_dot = u(4)*pi/180;% deg/s-->rad/s
    PosY = u(5);%单位为m
    PosX = u(6);%单位为米
    steer_L1 = u(7);
    steer_R1 = u(8);
    steer_deg = 0.5*(steer_L1 + steer_R1);
    delta_f_rad = steer_deg*pi/180;
    Beta = u(9)*pi/180;% 车辆质心侧偏角, Unit:deg-->rad
    slip_ratio_L1 = u(10);
    slip_ratio_L2 = u(11);
    slip_ratio_R1 = u(12);
    slip_ratio_R2 = u(13);
    Sf = 0.5*(slip_ratio_L1 +slip_ratio_R1);%前轮的滑移率
    Sr = 0.5*(slip_ratio_L2 +slip_ratio_R2);%后轮的滑移率
    alpha_L1 = u(14);
    alpha_L2 = u(15);
    alpha_R1 = u(16);
    alpha_R2 = u(17);
    alpha_f = 0.5*(alpha_L1 + alpha_R1);
    alpha_r = 0.5*(alpha_L2 + alpha_R2);
    Fy_l1      = round(10*u(18))/10; % N 
    Fy_l2      = round(10*u(19))/10; % N 
    Fy_r1      = round(10*u(20))/10; % N 
    Fy_r2      = round(10*u(21))/10; % N 
    Fyf        = Fy_l1 + Fy_r1;
    Fyr        = Fy_l2 + Fy_r2;
    Fx_L1      = u(22);
    Fx_L2      = u(23);
    Fx_R1      = u(24);
    Fx_R2      = u(25);    
    Fxf        = Fx_L1 + Fx_R1;
    Fxr        = Fx_L2 + Fx_R2;

    %-----Update augmented state vector--------------% 
    kesi    = zeros(Naug, 1);
    kesi(1) = y_dot;
    kesi(2) = x_dot;
    kesi(3) = PosPhi;
    kesi(4) = phi_dot;
    kesi(5) = PosY;
    kesi(6) = PosX;
    kesi(7) = delta_f_rad;   
    
    % C_cf C_cr C_lf C_lr分别为前后车轮的纵横向侧偏刚度，车辆固有参数    
    %----Estimate Lateral Cornering stiffness with RLS-------------------%  
    alpha_f_Hat = (Beta + phi_dot*VehiclePara.Lf/x_dot - delta_f_rad);
    [Fyf_hat, Ccf_1] = func_RLSEstimation_Ccf(alpha_f_Hat, Fyf);
    C_cf = sum(Ccf_1);
    if C_cf > -30000
        C_cf = -110000;
    end
    alpha_r_Hat = (Beta - phi_dot*VehiclePara.Lr/x_dot);
    [Fyr_hat, Ccr_1] = func_RLSEstimation_Ccr(alpha_r_Hat, Fyr);
    C_cr = sum(Ccr_1);
    if C_cr > -30000
        C_cr = -92000;
    end

    %-----Estimate Longitudinal Cornering stiffness with RLS--------------%
    Sf_Hat = Sf;
    [Fxf_hat, Clf_1] = func_RLSEstimation_Clf(Sf_Hat, Fxf);
    C_lf = sum(Clf_1);
    
    Sr_Hat = Sr;
    [Fxr_hat, Clr_1] = func_RLSEstimation_Clf(Sr_Hat, Fxr);
    C_lr = sum(Clr_1);
    
    %-----Use Constant tire stiffness  -------------------%  
%     C_cf = -57218; 
%     C_cr = -67587; 
%     
%     C_lf = 12650; 
%     C_lr = 99141;

    %*******Step(3): 更新状态方程 **********************************%    
    % 最基本也最重要的矩阵，是控制器的基础，采用动力学模型
    % 该矩阵与车辆参数密切相关，通过对动力学方程求解雅克比矩阵得到
    [Ad, Bd] = func_Model_linearization_Jacobian(kesi, Sf, Sr, C_cf, ...
                                             C_cr, C_lf, C_lr, ...
                                             MPCParameters, VehiclePara);

    A_cell = cell(2,2);
    A_cell{1,1} = Ad;
    A_cell{1,2} = Bd;
    A_cell{2,1} = zeros(Nu,Nx);
    A_cell{2,2} = eye(Nu);
    A = cell2mat(A_cell);
    
    B_cell = cell(2,1);
    B_cell{1,1} = Bd;
    B_cell{2,1} = eye(Nu);
    B = cell2mat(B_cell);

    C = [0 0 1 0 0 0 0;
         0 0 0 0 1 0 0];
                                                                             
	%*******Step(4): 参考轨迹生成 **********************************%    
%     %以下即为根据离散非线性模型预测下一时刻状态量 
%     [state_k1, Yita_ref] = func_Reftraj_doublelane(kesi, Sf, Sr, MPCParameters, ...
%                                          VehiclePara, C_cf, C_cr, C_lf, C_lr);
%     d_k = state_k1-Ad*kesi(1:6,1)-Bd*kesi(7,1);%根据falcone公式（2.11b）求得d(k,t)
%     d_piao_k = zeros(Nx+Nu, 1);%d_k的增广形式
%     d_piao_k(1:6,1) = d_k;
%     d_piao_k(7,1) = 0;
    
%      e_psi = kesi(3) - state_k1(3);
%      if(e_psi > pi)
%          e_psi = e_psi - 2*pi;
%      end
%      if(e_psi < -pi)
%          e_psi = e_psi + 2*pi;
%      end
%      e_y   = kesi(5) - state_k1(5);
     
     
    [WPIndex, Yita_ref, RefU] = func_RefTraj_LocalPlanning_DSL(MPCParameters,... 
                                            VehiclePara,... 
                                            WayPoints_IndexPre,... 
                                            Reftraj.DLS_path,... 
                                            kesi);
    if ( WPIndex <= 0)
       fprintf('Error: WPIndex <= 0 \n');    %出错
    else      
        WayPoints_IndexPre = WPIndex;        
    end
    
    d_piao_k = zeros(Nx+Nu, 1);%d_k的增广形式
   
    %****Step(5):  MPC formulation;********************% 
    %------Update prediction, ETA = PSI*kesi + GAMMA*PHI - Yref ----%
    PSI_cell=cell(Np,1);
    THETA_cell=cell(Np,Nc);
    GAMMA_cell=cell(Np,Np);
    PHI_cell=cell(Np,1);
    for p=1:1:Np;
        PHI_cell{p,1}=d_piao_k;%理论上来说，这个是要实时更新的，但是为了简便，这里又一次近似
        for q=1:1:Np;
            if q<=p;
                GAMMA_cell{p,q}=C*A^(p-q);
            else 
                GAMMA_cell{p,q}=zeros(Ny,Nx+Nu);
            end 
        end
    end
    for j=1:1:Np
     PSI_cell{j,1}=C*A^j;
        for k=1:1:Nc
            if k<=j
                THETA_cell{j,k}=C*A^(j-k)*B;
            else 
                THETA_cell{j,k}=zeros(Ny,Nu);
            end
        end
    end
    PSI=cell2mat(PSI_cell);%size(PSI)=[Ny*Np Nx+Nu]
    THETA=cell2mat(THETA_cell);%size(THETA)=[Ny*Np Nu*Nc]
    GAMMA=cell2mat(GAMMA_cell);%大写的GAMMA
    PHI=cell2mat(PHI_cell);
    
    %------Update Q and R
    temp = [CostWeights.Wephi, CostWeights.Wey];
    Qq = diag(temp);
    Q = kron(eye(Np), Qq);
    R = kron(eye(Nc), CostWeights.WDdeltaf);
    
    %------Update H and f, J=0.5*DU'*H*DU + f'*DU
    H = THETA'*Q*THETA + R;
    H = 0.5*(H+H');
    error_1 = PSI*kesi + GAMMA*PHI - Yita_ref; %求偏差
    f = error_1'*Q*THETA;
    g = f';
    
    %------Update Constaints and bounds 生成约束---%
    %控制量约束
    A_t = zeros(Nc,Nc);%见falcone论文 P181
    for p = 1:1:Nc
        for q = 1:1:Nc
            if (p >= q)
                A_t(p,q) = 1;
            else 
                A_t(p,q) = 0;
            end
        end 
    end 
    A_I = kron(A_t,eye(Nu));%求克罗内克积
    
    Ut=kron(ones(Nc,1), delta_f_rad);
    umin = -Constraints.umax;%维数与控制变量的个数相同
    umax = Constraints.umax; 
    Umin=kron(ones(Nc,1), umin);
    Umax=kron(ones(Nc,1), umax);

    A_cons_cell = { A_I; 
                    -A_I};
    A_cons=cell2mat(A_cons_cell);%（求解方程）状态量不等式约束增益矩阵，转换为绝对值的取值范围  
   
    b_cons_cell = { Umax - Ut; 
                    -Umin + Ut};
    b_cons = cell2mat(b_cons_cell);%（求解方程）状态量不等式约束的取值
  
%     ycmax = Constraints.ycmax;
%     ycmin = Constraints.ycmin;
%     Ycmax = kron(ones(Np,1),ycmax);
%     Ycmin = kron(ones(Np,1),ycmin);  %输出量约束   
%    
%     A_cons_cell = { A_I; 
%                     -A_I;
%                     THETA;
%                     -THETA};
%     A_cons=cell2mat(A_cons_cell);%（求解方程）状态量不等式约束增益矩阵，转换为绝对值的取值范围  
%    
%     b_cons_cell = { Umax - Ut; 
%                     -Umin + Ut; 
%                     Ycmax - error_1; 
%                     -Ycmin + error_1};
%     b_cons = cell2mat(b_cons_cell);%（求解方程）状态量不等式约束的取值
   
    lb=kron(ones(Nc,1), -Constraints.dumax);
    ub=kron(ones(Nc,1), Constraints.dumax);
    
    %****Step(9):  Call quadprog for MPC solver;********************% 
%     options = optimset('quadprog', 'Algorithm', 'active-set');
    options = optimoptions('quadprog', 'Display','off', ...
                            'Algorithm', 'active-set'); 
    DU0 = WarmStart;
    [DU, FVAL, EXITFLAG] = quadprog(H, g, A_cons, b_cons, [], [], lb, ub, DU0, options); %
%     [DU, FVAL, EXITFLAG] = quadprog(H, g, [], [], [], [], lb, ub, DU0, options); %
    WarmStart = shiftHorizon(DU, Nu);     % Prepare restart, nominal close loop 
    
    Steer_SW_deg = 18 * (delta_f_rad + DU(1))*180/pi;
    
    t_Elapsed = toc( t_Start ); %computation time
end % end of if Initialflag < 2 % 

sys = [Steer_SW_deg; t_Elapsed; PosX; PosY; PosPhi; e_psi; e_y];
% end  %End of mdlOutputs.

%==============================================================
% sub functions
%==============================================================  
function U0 = shiftHorizon(U, Nu) %shift control horizon
    U0 = [U(Nu+1:size(U,1)); zeros(Nu,1)]; % shiftHorizon：Prepare restart
    
function [WPIndex, Yita_ref, RefU] = func_RefTraj_LocalPlanning_DSL( MPCParameters, VehiclePara, WayPoints_Index, WayPoints_Collect, VehStateVector)
    lf = VehiclePara.Lf;
    lr = VehiclePara.Lr;
    lfr = VehiclePara.L;
    m   = VehiclePara.m;
    Iz  = VehiclePara.Iz;   %I为车辆绕Z轴的转动惯量，车辆固有参数  
    Ts = MPCParameters.Ts;
    Np  = MPCParameters.Np;
    Nx  = MPCParameters.Nx;
    Nu  = MPCParameters.Nu;
    
    PosPsi  = VehStateVector(3);  
    PosY    = VehStateVector(5);
    PosX    = VehStateVector(6);
%*********** WaypointData2VehicleCoords ************************% 
    ds          = 0.1;%m
    WPNum       = length(WayPoints_Collect(:,1));
    
    %--------先找到参考路径上距离车辆最近的点--------------------------%  
    Dist_MIN    = 1000;
    index_min   = 0;
    for i=WayPoints_Index:1:WPNum 
        deltax  = WayPoints_Collect(i,1) - PosX;
        deltay  = WayPoints_Collect(i,2) - PosY;
        Dist    = sqrt(power(deltax,2) + power(deltay,2)); %路点到车辆重心的距离
        if Dist < Dist_MIN
            Dist_MIN = Dist; 
            index_min = i;
        end
    end
    
    if (index_min < 1) 
        WPIndex = -1; %如果没有找到则。。
    else if ( index_min >= WPNum)
            WPIndex = -2; %如果没有找到则。。
        else
            WPIndex = index_min;
        end
    end

    Yita_ref_cell=cell(Np,1);
    for p=1:1:Np
%         X_ref = WayPoints_Collect(WPIndex+p, 1);
        Y_ref = WayPoints_Collect(WPIndex+p, 2);
        Heading_ref = WayPoints_Collect(WPIndex+p, 2);
        Yita_ref_cell{p,1} = [Heading_ref; Y_ref];
    end
    Yita_ref=cell2mat(Yita_ref_cell);
    
    RefU = zeros(Np,1);
    
% end % End of func

function [state_k1, Yita_ref] = func_Reftraj_doublelane(kesi, Sf, Sr, MPCParameters, VehiclePara, Ccf, Ccr, Clf, Clr)

% 双移线轨迹形状参数
    shape=2.4;%参数名称，用于参考轨迹生成
    dx1=25; dx2=21.95;%没有任何实际意义，只是参数名称
    dy1=4.05; dy2=5.7;%没有任何实际意义，只是参数名称
    Xs1=27.19; Xs2=56.46;%参数名称

    lf = VehiclePara.Lf;
    lr = VehiclePara.Lr;
    lfr = VehiclePara.L;
    m   = VehiclePara.m;
    Iz  = VehiclePara.Iz;   %I为车辆绕Z轴的转动惯量，车辆固有参数  
    Ts = MPCParameters.Ts;
    Np  = MPCParameters.Np;
    Nx  = MPCParameters.Nx;
    Nu  = MPCParameters.Nu;
   
    y_dot   = kesi(1);%u(1)==X(1)
    x_dot   = kesi(2);%u(2)==X(2)
    phi     = kesi(3); %u(3)==X(3)
    phi_dot = kesi(4);
    Y       = kesi(5);
    X       = kesi(6);
    delta_f = kesi(7);   

    state_k1 = zeros(Nx, 1);
    state_k1(1,1)=y_dot+Ts*(-x_dot*phi_dot+2*(Ccf*(delta_f-(y_dot+lf*phi_dot)/x_dot)+Ccr*(lr*phi_dot-y_dot)/x_dot)/m);
    state_k1(2,1)=x_dot+Ts*(y_dot*phi_dot+2*(Clf*Sf+Clr*Sr+Ccf*delta_f*(delta_f-(y_dot+phi_dot*lf)/x_dot))/m);
    state_k1(3,1)=phi+Ts*phi_dot;
    state_k1(4,1)=phi_dot+Ts*((2*lf*Ccf*(delta_f-(y_dot+lf*phi_dot)/x_dot)-2*lr*Ccr*(lr*phi_dot-y_dot)/x_dot)/Iz);
    state_k1(5,1)=Y+Ts*(x_dot*sin(phi)+y_dot*cos(phi));
    state_k1(6,1)=X+Ts*(x_dot*cos(phi)-y_dot*sin(phi));
    
%     T_all = 20.0;    
    Yita_ref_cell=cell(Np,1);
    X_DOT=x_dot*cos(phi)-y_dot*sin(phi);%惯性坐标系下纵向速度
    for p=1:1:Np
        X_predict = X+X_DOT*p*Ts;%首先计算出未来X的位置
        z1      = shape/dx1*(X_predict-Xs1) - shape/2;
        z2      = shape/dx2*(X_predict-Xs2) - shape/2;
        Y_ref   = dy1/2*(1+tanh(z1)) - dy2/2*(1+tanh(z2));
        phi_ref = atan(dy1*(1/cosh(z1))^2*(1.2/dx1) - dy2*(1/cosh(z2))^2*(1.2/dx2));
        Yita_ref_cell{p,1} = [phi_ref; Y_ref];
    end
    Yita_ref=cell2mat(Yita_ref_cell);
% end % End of func















