function [sys,x0,str,ts] =Main_RollDynamicsEstimation_UKF_Modified(t,x,u,flag)
%***************************************************************%
% This dynamics state estimation is based on "Roll Prediction-based Optimal
% Control for Safe Path Following" by Sanghyun Hong and J. Karl Hedrick.
%  对于Vy的估计效果不好，可能是因为动力学模型不够精确
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
sizes.NumDiscStates  = 5;  %模块离散状态变量的个数,实际上没有用到这个数值，只是用这个来表示离散模块
sizes.NumOutputs     = 10;  %S函数的输出，包括控制量和其它监测量
sizes.NumInputs      = 39; %S函数模块输入变量的个数，即CarSim的输出量
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

global Initial_State; 
global P; 
global Numx;
Initial_State=[0; 0 ; 0 ; 0; 0];  % initial state
Numx = numel(Initial_State);
P = eye(Numx);                 % initial state covraiance

% global PreviousYawrate;
% PreviousYawrate = 0;

global Wm;
global Wc;
global r_sigmas;
[Wm, Wc, r_sigmas] = UTSetup(Numx);

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
global Initial_State; 
global P; 
% global PreviousYawrate;
global Wm;
global Wc;
global r_sigmas;
global Numx;

Ts  = 0.05;
Iz  = 2059.2;
m   = 1600;
ms  = 1430;
g   = 9.81;
hs  = 0.65;
Lf  = 1.12;  %a
Lr  = 1.48;  %b
% L   = 2.6;
Lw  = 1.565;
Ix  = 700.7;

rho     = 1.206;
Cd      = 0.15;
Af      = 3.1;
f_roll  = 0.02;
Kroll   = 145330; %42075;
Croll   = 4500.5; %5737.5;

%% Updatd measured State of Vehicle
    [VehStateMeasured, ParaHAT] = func_StateEstimation(u); % unpack data from Carsim   

    Vx          = VehStateMeasured.x_dot; 
    Vy          = VehStateMeasured.y_dot; 
    Yawrate     = VehStateMeasured.phi_dot; % rad/s
    Ax          = VehStateMeasured.Ax; % x_dot
    Ay          = VehStateMeasured.Ay; % y_dot
    fwa         = VehStateMeasured.fwa;
    Beta        = VehStateMeasured.beta;%rad    
    if 0 == Vx
        Vx = 10;
    end
    AyG_SM      = ParaHAT.AyG_SM;    
    Ay_Bf_SM    = ParaHAT.Ay_Bf_SM;    
    Ax_SM       = ParaHAT.Ax_SM; %m/s^2
    alphaf      = ParaHAT.alphaf;
    alphar      = ParaHAT.alphar;
    
    Rollangle   = ParaHAT.Roll;
    Rollrate    = ParaHAT.Rollrate;
    
    Faero   = 0.5*rho*Cd*Af*Vx^2;
    %     Fyf         = ParaHAT.Fyf;
    %     Fyr         = ParaHAT.Fyr;   

    Fcfl = ParaHAT.Fy_l1;
    Fcfr = ParaHAT.Fy_r1;
    Fcrl = ParaHAT.Fy_l2;
    Fcrr = ParaHAT.Fy_r2;

    Flfl = ParaHAT.Fx_L1;
    Flfr = ParaHAT.Fx_R1;
    Flrl = ParaHAT.Fx_L2;
    Flrr = ParaHAT.Fx_R2;

    Fxfl = Flfl*cos(fwa) -Fcfl*sin(fwa);
    Fxfr = Flfr*cos(fwa) -Fcfr*sin(fwa);
    Fxrl = Flrl;
    Fxrr = Flrr;
    Fyfl = Flfl*sin(fwa) + Fcfl*cos(fwa);
    Fyfr = Flfr*sin(fwa) + Fcfr*cos(fwa);
    Fyrl = Fcrl;
    Fyrr = Fcrr;
    Fx   = Fxfl+Fxfr+Fxrl+Fxrr;
    Fy   = Fyfl+Fyfr+Fyrl+Fyrr;
    


%% Dual UKF process for estimation
%-----------------Vehicle Dynamics Model according to "Sanghyun Hong and J. Karl Hedrick"-------%
% Let s = [Vx; Vy; Yawrate; Rollangle; Rollrate];
%     x1_dot = (Fx-Faero)/m + Vy*Yawrate;
%     x2_dot = Fy/m - Vx*Yawrate;    
%     x3_dot = ( Lf*(Fyfl+Fyfr)-Lr*(Fyrl+Fyrr) + Lw*(Fxfr+Fxrr-Fxfl-Fxrl)/2 )/Iz;
%     x4_dot = rollrate;
%     x5_dot = ( ms*hs*( Fy/m+ g*Rollangle ) - Kroll*Rollangle - Croll*Rollrate)/Ix;
% Discrete the previous model with Ts, then we get:
fstate = @(s)[  s(1) + Ts*((Fx-Faero-f_roll*m*g)/m + s(2)*s(3)) - hs*s(3)*s(5);  %  2*
                s(2) + Ts*(Fy/m - s(1)*s(3)) + hs* (( m*hs*( Fy/m + g*s(4)) -Kroll*s(4) - Croll*s(5) )/Ix); % 
                s(3) + Ts*(( Lf*(Fyfl+Fyfr) - Lr*(Fyrl+Fyrr) + Lw*(Fxfr+Fxrr-Fxfl-Fxrl)/2 )/Iz);
                s(4) + Ts*s(5);
                s(5) + Ts*(( m*hs*( Fy/m + g*s(4)) -Kroll*s(4) - Croll*s(5) )/Ix) ]; 
%------------- ay_sensor, Ax, yawrate and Vx can be measured -----------*%
Ymeasurement = [AyG_SM; Vx; Yawrate; Ax_SM];  
Ny = 4; 
hmeas = @(s)[ Fy/m  + g*s(4); % - hs* (( m*hs*( Fy/m + g*s(4)) -Kroll*s(4) - Croll*s(5) )/Ix)
              s(1);
              s(3);
             (Fx-Faero)/m - f_roll*g]; 
         
%-------------------------Dynamics of DualUKF----------------------------% 
Q = 0.1 * eye(Numx);
R = 0.1 * eye(Ny);
 
X=Sigmas(Initial_State,P,r_sigmas);                           %sigma points around x
[x1,X1,P1,X2]=UnscentedTransform(fstate, X, Wm, Wc, Numx, Q);      %unscented transformation of process

[z1,Z1,P2,Z2]=UnscentedTransform(hmeas, X1, Wm, Wc, Ny, R);       %unscented transformation of measurments

    % Update after the measurement of 
     P12=X2*diag(Wc)*Z2';                        %transformed cross-covariance   
%     %********************Normal solution********************%
%     K=P12*inv(P2);
%     X0=x1+K*(Ymeasurement-z1);                              %state update
%     P0=P1-K*P12';                                %covariance update
     %********************Cholesky solution********************%
    R=chol(P2);            %Cholesky factorization
    U=P12/R;                    %K=U/R'; Faster because of back substitution
    X0 = x1+U*(R'\(Ymeasurement-z1));     %Back substitution to get state update
    P0=P1-U*U';                   %Covariance update, U*U'=P12/R/R'*P12'=K*P12. 
        
Initial_State = X0;
P             = P0;
% PreviousYawrate   = Yawrate;

sys = [X0; Vx; Vy; Yawrate; Rollangle; Rollrate];

% end  %End of mdlOutputs.

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
    VehStatemeasured.delta_f = (round(10*0.5*(ModelInput(8)+ ModelInput(9)))/10)*pi/180; % deg-->rad
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
    
    HATParameter.GearStat    = ModelInput(30);
    HATParameter.Roll        = ModelInput(31)*pi/180;% deg-->rad 
    HATParameter.Rollrate    = ModelInput(32)*pi/180;% deg/s-->rad/s
    HATParameter.Roll_accel  = ModelInput(33); % rad/s^2
    HATParameter.Z0          = ModelInput(34); %m
%     VehStatemeasured.Station     = ModelInput(35); %m
    HATParameter.Zcg_TM      = ModelInput(35); %m
    HATParameter.Zcg_SM      = ModelInput(36); %m
    HATParameter.AyG_SM       = ModelInput(37)*g; %m/s^2, acceleration measured by accelerometer
    HATParameter.Ay_Bf_SM    = ModelInput(38)*g; %m/s^2
    HATParameter.Ax_SM       = ModelInput(39)*g; %m/s^2
   
function [Wm, Wc, r_sigmas] = UTSetup(L)
    alpha=1e-2;                                 %default, tunable
    ki=1;                                       %default, tunable, generally set to 3-L
    beta=2;                                     %default, tunable
    lambda=alpha^2*(L+ki)-L;                    %scaling factor
    c=L+lambda;                                 %scaling factor
    Wm=[lambda/c 0.5/c+zeros(1,2*L)];           %weights for means
    %Wm=[lambda/c ones(1,2*L)/(2*c)];           %weights for means
    Wc=Wm;                                      %length=2*L+1
    Wc(1)=Wc(1)+(1-alpha^2+beta);               %weights for covariance
    r_sigmas=sqrt(c);


function [y,Y,P,Y1] = UnscentedTransform(f,X,Wm,Wc,n,R)
%Unscented Transformation
%Input:
%        f: nonlinear map
%        X: sigma points
%       Wm: weights for mean
%       Wc: weights for covraiance
%        n: numer of outputs of f
%        R: additive covariance
%Output:
%        y: transformed mean
%        Y: transformed smapling points
%        P: transformed covariance
%       Y1: transformed deviations

L=size(X,2);
y=zeros(n,1);
Y=zeros(n,L);
for k=1:L                   
    Y(:,k)=f(X(:,k));       
    y=y+Wm(k)*Y(:,k);       
end
Y1=Y-y(:,ones(1,L));
P=Y1*diag(Wc)*Y1'+R;          

function X = Sigmas(x,P,c)
%Sigma points around reference point
%Inputs:
%       x: reference point
%       P: covariance
%       c: coefficient
%Output:
%       X: Sigma points

A = c*chol(P)';
Y = x(:,ones(1,numel(x)));
X = [x Y+A Y-A]; 



 