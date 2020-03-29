function    [Ad, Bd] = func_Model_linearization_Jacobian(kesi, Sf, Sr, Ccf, Ccr, Clf, Clr, MPCParameters, VehiclePara)
%***************************************************************%
% 根据简化动力学模型(考虑小角度假设下)计算雅克比矩阵
% 计算出的雅克比矩阵是与车辆固有参数紧密相关
%---------------------------------------------------------------%
% Published by: Kai Liu
% Email:leoking1025@gmail.com
% My homepage: https://sites.google.com/site/kailiumiracle/  
%***************************************************************%
    %----------车辆参数定义 -----------%
    syms y_dot x_dot phi phi_dot Y X;%车辆状态量
    syms delta_f  %前轮偏角,控制量
    
    Ts = MPCParameters.Ts;
    a  = VehiclePara.Lf;
    b  = VehiclePara.Lr;
    m  = VehiclePara.m;
    Iz = VehiclePara.Iz;

    %----车辆动力学模型-------------%
    dy_dot = -x_dot*phi_dot + 2*(Ccf*((y_dot+a*phi_dot)/x_dot - delta_f) + Ccr*(y_dot - b*phi_dot)/x_dot)/m;
    dx_dot = y_dot*phi_dot + 2*(Clf*Sf + Clr*Sr + Ccf*((y_dot + phi_dot*a)/x_dot - delta_f)*delta_f)/m;
    dphi_dot = (2*a*Ccf*((y_dot+a*phi_dot)/x_dot - delta_f) - 2*b*Ccr*(y_dot - b*phi_dot)/x_dot)/Iz;
    Y_dot = x_dot*sin(phi) + y_dot*cos(phi);
    X_dot = x_dot*cos(phi) - y_dot*sin(phi);

    %----雅克比矩阵求解-------------%
    Dynamics_func = [dy_dot; dx_dot; phi_dot; dphi_dot; Y_dot; X_dot];%动力学模型
    state_vector = [y_dot,x_dot,phi,phi_dot,Y,X];%系统状态量
    control_input = delta_f;
    A_t = jacobian(Dynamics_func, state_vector);  %矩阵A(t)-连续
    B_t = jacobian(Dynamics_func, control_input); %矩阵B(t)-连续

    %----将连续矩阵转换为离散矩阵-------------%
    % 采用Forward Euler Method近似算法  A = I+Ts*A(t),B = Ts*B(t)
    I_6 = eye(6);
    Ad_temp = I_6 + Ts * A_t;
    Bd_temp = Ts * B_t;

    %----提取车辆状态参数-------------%
    y_dot   = kesi(1);
    x_dot   = kesi(2);
    phi     = kesi(3);
    phi_dot = kesi(4);
    Y       = kesi(5);
    X       = kesi(6);
    delta_f = kesi(7);
    
    Ad = eval(Ad_temp);
    Bd = eval(Bd_temp);
    
end % end of func.