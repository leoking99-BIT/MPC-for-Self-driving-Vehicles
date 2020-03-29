function [sys,x0,str,ts] = Main_MPC_GivenPathTracking(t,x,u,flag)
%***************************************************************%
%   基于车辆运动学模型实现给定参考轨迹的跟踪
%   控制量为速度和前轮偏角
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
% Initialization
%==============================================================
function [sys,x0,str,ts] = mdlInitializeSizes
%***************************************************************%
% Call simsizes for a sizes structure, fill it in, and convert it 
% to a sizes array.
%***************************************************************%
sizes = simsizes;
sizes.NumContStates  = 0;
sizes.NumDiscStates  = 3;
sizes.NumOutputs     = 9;
sizes.NumInputs      = 10;
sizes.DirFeedthrough = 1; % Matrix D is non-empty.
sizes.NumSampleTimes = 1;
sys = simsizes(sizes); 
x0 = zeros(sizes.NumDiscStates,1);%initial the  state vector
str = [];             % Set str to an empty matrix.
ts  = [0.05 0];       % sample time: [period, offset]

global MPCParameters; 
    MPCParameters.Ts      = 0.05; %Set the sample time
    MPCParameters.Np      = 30;% predictive horizon
    MPCParameters.Nc      = 30;% control horizon
    MPCParameters.Nx      = 3; %number of state variables
    MPCParameters.Nu      = 2; %number of control inputs
    MPCParameters.Ny      = 3; %number of output variables  
    MPCParameters.Qx      = 10; % cost weight factor 
    MPCParameters.Qy      = 10; % cost weight factor 
    MPCParameters.Qphi    = 10; % cost weight factor 
    MPCParameters.R_dv    = 5; % cost weight factor 
    MPCParameters.R_du    = 5; % cost weight factor 
    MPCParameters.qp_solver = 0; %0: default, quadprog; 1:qpOASES; 2:CVXGEN
    MPCParameters.refspeedT = 0; %0: default, straight-lane profile; 
                                 %1:circle path profile
    MPCParameters.umin      = [-0.2; -0.436;];  % 
    MPCParameters.umax      = [0.2;  0.436];  % the max 
    MPCParameters.dumin     = [-0.05; -0.0082;]; % 
    MPCParameters.dumax     = [0.05; 0.0082]; % 

global WarmStart;
    WarmStart = zeros(MPCParameters.Nu * MPCParameters.Nc,1);
global InitialGapflag; 
    InitialGapflag = 0; % Ignore the first few inputs from CarSim
% Initialize the discrete states.

%End of mdlInitializeSizes
		      
%==============================================================
% Update the discrete states
%==============================================================
function sys = mdlUpdates(t,x,u)
%  目前没有用到这个过程；在后期的程序模块化时可以继续开发这个功能。
    sys = x;
% end     %End of mdlUpdate.

%==============================================================
% Calculate outputs
%==============================================================
function sys = mdlOutputs(t,x,u)
    
global InitialGapflag;
global MPCParameters;
global WarmStart;

Ts = MPCParameters.Ts;
Np = MPCParameters.Np;
Nc = MPCParameters.Nc;
Nx = MPCParameters.Nx;
Nu = MPCParameters.Nu;

% lfr = 2.78; % D-class Sedan
lfr = 2.6; % D-class SUV

Vel_des   =  0.0;
Steer_des =  0.0;
            
t_Start = tic; % 开始计时 
if InitialGapflag < 2 %  get rid of the first two inputs
    InitialGapflag = InitialGapflag + 1;%
else
    %***********Step (1). Update vehicle states *************************% 
    % 提取CarSim输入到Simulink的数据
    x_L2 = u(1); %左后轮x坐标
    x_R2 = u(2); %右后轮x坐标
    y_L2 = u(3); %左后轮y坐标
    y_R2 = u(4); %右后轮y坐标   
    Yaw_angle  = u(5)*pi/180;%航向角Unit：deg-->rad
    Steer_SW = u(6); %方向盘角度
    Steer_L1 = u(7); %左前轮偏角
    Steer_R1 = u(8); %右前轮偏角
    Vx_L2 = u(9);  %左后轮纵向速度，Unit:km/h
    Vx_R2 = u(10); %右后轮纵向速度，Unit:km/h
    
    VehPos.X = 0.5*(x_L2 + x_R2);%后轴中心X坐标，Unit:m
    VehPos.Y = 0.5*(y_L2 + y_R2);%后轴中心Y坐标，Unit:m
    VehPos.Yaw_angle = Yaw_angle;%车辆航向角，Unit:rad
    
    Vx_km_h = 0.5*(Vx_L2 + Vx_R2);%后轴中心处纵向速度,Unit：km/h
    Steer_deg = 0.5*(Steer_L1 + Steer_R1);%等效前轮偏角，Unit：deg
    Vx_m_s  = Vx_km_h/3.6;%%后轴中心处纵向速度 in (m/s),Unit：m/s    
    Steer_rad = Steer_deg*pi/180;%等效前轮偏角in (rad)，Unit：degs-->rad;
    
    %********Step(2): Generate reference speed profile *******************%
    switch MPCParameters.refspeedT,
        case 0 % default, straight lane profile
            %----设定直道参考道路的起点和终点----------------------%
            V_ref = 5; %参考车速， unit:m/s
            StartPoint = [-1, 5]; %起点的x坐标和y坐标
            EndPoint = [10, 5]; %起点的x坐标和y坐标
            PrjPoint = func_GetProjectPoint_StraightLane(VehPos, StartPoint, EndPoint);
        case 1 % circle lane profile
            %----设置Circle形式 的参考道路曲线--------------%
            V_ref = 5; %参考车速， unit:m/s
            Circle_center = [0, 25]; %圆心的横纵向坐标, unit:m
            Circle_radius = 5; %参考车速，unit:m/s
            PrjPoint = func_GetProjectPoint_Circle(VehPos, lfr, Circle_center, Circle_radius);
        otherwise % Unexpected flags %
            error(['unexpected path-profile:',num2str(MPCParameters.refspeedT)]); % Error handling
    end %  end of switch 
    
    %****Step(3): update lateral vehilce model***%
    kesi=zeros(Nx + Nu, 1);
    kesi(1)=VehPos.X - PrjPoint.X;%
    kesi(2)=VehPos.Y- PrjPoint.Y;%
    kesi(3)=VehPos.Yaw_angle - PrjPoint.Yaw_angle; %
    kesi(4)=Vx_m_s;
    kesi(5)=Steer_rad;
    Um = [Vx_m_s; Steer_rad];
   
    cos2 = cos(PrjPoint.steer_rad)^2;
    Ad = [ 1    0   -V_ref*sin(PrjPoint.Yaw_angle)*Ts;
           0    1   V_ref*cos(PrjPoint.Yaw_angle)*Ts;
           0    0   1;];
    Bd = [ cos(PrjPoint.Yaw_angle)*Ts           0;
           sin(PrjPoint.Yaw_angle)*Ts           0;
           tan(PrjPoint.steer_rad)*Ts/lfr       V_ref*Ts/(lfr*cos2);];
    
    A_aug_cell = cell(2,2);
    A_aug_cell{1,1} = Ad;
    A_aug_cell{1,2} = Bd;
    A_aug_cell{2,1} = zeros(Nu, Nx);
    A_aug_cell{2,2} = eye(Nu);
    A = cell2mat(A_aug_cell);
 
    B_aug_cell = cell(2,1);
    B_aug_cell{1,1} = Bd;
    B_aug_cell{2,1} = eye(Nu);
    B = cell2mat(B_aug_cell);
    
    C=[ 1 0 0 0 0;
        0 1 0 0 0;
        0 0 1 0 0];
    
    %****Step(4):  MPC formulation;********************%
    %Update Theta and PHI for future states prediction
    PHI_cell=cell(Np,1);
    THETA_cell=cell(Np,Nc);
    for j=1:1:Np
        PHI_cell{j,1}=C*A^j; %Ny * N_aug=Nx+Nu
        for k=1:1:Nc
            if k <= j
                THETA_cell{j,k}=C*A^(j-k)*B; %Ny * Nu
            else 
                THETA_cell{j,k}=zeros(Nx,Nu);%
            end
        end
    end
    PHI=cell2mat(PHI_cell);%size(PHI)=[Nx*Np, Nx+Nu]
    THETA=cell2mat(THETA_cell);%size(THETA)=[Nx*Np, Nu*Nc)]    
    
    %Update H and f for cost function J = 0.5U'HU +　f'U
    Qq = diag([MPCParameters.Qx, MPCParameters.Qy, MPCParameters.Qphi]);
    Rr = diag([MPCParameters.R_dv, MPCParameters.R_du]);
    Q = kron(eye(Np), Qq); % Nx*Np,Nx*Np
    R = kron(eye(Nc), Rr); % Nu*Nc, Nu*Nc
    
    H = THETA'*Q*THETA + R;

    E_error = PHI*kesi;
    g       = THETA' * Q * E_error;
    
    %Update bounds and constraints 
    [A, b, Aeq, beq, lb, ub] = func_Constraints_du_quadprog(MPCParameters, Um);
        
    %****Step(5):  Call qp-solver********************%
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
    WarmStart = shiftHorizon(U, Nu);     % Prepare restart, nominal close loop 
    if (1 ~= EXITFLAG) %if optimization NOT succeeded.
        U(1) = 0.0;
        U(2) = 0.0;
        fprintf('MPC solver not converged!\n');                  
    end
    Vel_des   =  U(1) +  kesi(4) + V_ref;
    Steer_des = U(2) +  kesi(5) + PrjPoint.steer_rad;
    Steer_des =  18.0 * (Steer_des*180/pi);

    
end % End of if Initialflag < 1 %

t_Elapsed = toc( t_Start ); %computation time 

sys= [Vel_des; Steer_des; t_Elapsed; Car_X; Car_Y; Yaw; Vx_m_s; Steer_rad; Steer_SW];
% end % End of mdlOutputs

%==============================================================
% sub functions
%============================================================== 
function [PPx,PPy,slope]=func_GetProjectPoint(Xb,Yb,Xn,Yn,Xc,Yc)
    DifX=Xn-Xb;
    DifY=Yn-Yb;
    if DifX == 0
        x=Xn;
        y=Yc;
    else if DifY == 0
            x=Xc;
            y=Yn;
        else  %k!=0
            Kindex = DifY/DifX;
            bindex = Yb-Kindex*Xb;
            K = -1/Kindex;
            b = Yc-K*Xc;
            x = (bindex-b)/(K-Kindex);
            y = K*x+b;
        end     
    end
    de=atan2(DifY,DifX);
    if de<0
        de=de+2*pi;
    end   
    PPx=x;
    PPy=y;
    slope = de;
% end % End of func

function PrjPoint = func_GetProjectPoint_StraightLane(VehPos, StartPoint, EndPoint)

   [PrjPoint.X,PrjPoint.Y,PrjPoint.Yaw_angle] = func_GetProjectPoint(StartPoint(1), StartPoint(2), ...
                                          EndPoint(1), EndPoint(2), ...
                                          VehPos.X, VehPos.Y);
   PrjPoint.steer_rad = 0.0;
% end % End of func

function PrjPoint = func_GetProjectPoint_Circle(VehPos, lfr, Circle_center, Circle_radius)
    CEN_x = Circle_center(1);
    CEN_y = Circle_center(2);

    Alpha_init = Func_Alpha_Pos(CEN_x,CEN_y,VehPos.X, VehPos.Y);%首先根据车辆位置和圆心确定alpha
    PrjPoint.X = Circle_radius*cos(Alpha_init) + CEN_x;%x
    PrjPoint.X = Circle_radius*sin(Alpha_init) + CEN_y;%y
    PrjPoint.Yaw_angle = Func_Theta_Pos(Alpha_init);%theta  
    
    PrjPoint.steer_rad = atan(lfr/Circle_radius);
%end %ENd of func
    
function U0 = shiftHorizon(U, Nu) %shift control horizon
    U0 = [U(Nu+1:size(U,1)); zeros(Nu,1)]; % shiftHorizon：Prepare restart
    
function  [A, b, Aeq, beq, lb, ub] = func_Constraints_du_quadprog(MPCParameters, um)
%************************************************************************%
% generate the constraints of the vehicle
%  
%************************************************************************%
    Np   = MPCParameters.Np;
    Nc   = Np;
    Nu   = MPCParameters.Nu;
    dumin = MPCParameters.dumin;
    dumax = MPCParameters.dumax;
    umin = MPCParameters.umin;  
    umax = MPCParameters.umax;  
    Umin = kron(ones(Nc,1),umin);
    Umax = kron(ones(Nc,1),umax);
    Ut   = kron(ones(Nc,1),um);
%----(1) A*x<=b----------%
    A_t=zeros(Nc,Nc);
    for p=1:1:Nc
        for q=1:1:Nc
            if p >= q 
                A_t(p,q)=1;
            else 
                A_t(p,q)=0;
            end
        end 
    end
    A_I=kron(A_t,eye(Nu));%对应于falcone论文约束处理的矩阵A,求克罗内克积
    
    A_cell=cell(2,1);
    A_cell{1,1} = A_I; %
    A_cell{2,1} = -A_I;
    A=cell2mat(A_cell);  %
    
    b_cell=cell(2, 1);
    b_cell{1,1} = Umax - Ut; %
    b_cell{2,1} = -Umin + Ut;
    b=cell2mat(b_cell);  % 

%----(2) Aeq*x=beq----------%
    Aeq = [];
    beq = [];

%----(3) lb=<x<=ub----------%
    lb=kron(ones(Nc,1), dumin);
    ub=kron(ones(Nc,1), dumax);
% end %EoF

function [A_I, lb, ub, lbA, ubA] = func_Constraints_du_qpOASES(MPCParameters, um)
    Np   = MPCParameters.Np;
    Nc   = Np;
    dumin = MPCParameters.dumin;
    dumax = MPCParameters.dumax;
    umin = MPCParameters.umin;
    umax = MPCParameters.umax;  
    Umin = kron(ones(Nc,1), umin);
    Umax = kron(ones(Nc,1), umax);
    Ut   = kron(ones(Nc,1),um);
%----(1) lbA <= A_t*x<=ubA----------%
    A_t=zeros(Nc,Nc);
    for p=1:1:Nc
        for q=1:1:Nc
            if p >= q 
                A_t(p,q)=1;
            else 
                A_t(p,q)=0;
            end
        end 
    end 
    A_I=kron(A_t,eye(Nu));%对应于falcone论文约束处理的矩阵A,求克罗内克积
    ubA = Umax - Ut; %
    lbA = Umin - Ut;
%---- lb=<x<=ub----------%
    lb=kron(ones(Nc,1), dumin);
    ub=kron(ones(Nc,1), dumax);
% end %EoF





