function [Updated_state] = func_UpdateState_dsolve_2_7(Previous_States, lfr, Vx_m_s, Steer_rad, Ts)
% validation of dynamic model (2.26) use dsolve
    
    X_init      = Previous_States.X_pred;
    Y_init      = Previous_States.Y_pred;
    Yaw_init    = Previous_States.Yaw_pred;

    XOUT=dsolve('Dx-Vx_m_s*cos(z) = 0', ...
                'Dy-Vx_m_s*sin(z) = 0', ...
                'Dz-Vx_m_s*tan(Steer_rad)/lfr = 0', ...
                'x(0) = X_init', ...
                'y(0) = Y_init', ...
                'z(0) = Yaw_init');%z denotes yawangle
    t=0.05; % independent variable is 't' means sample time
    X_pred   = eval(XOUT.x);
    Y_pred   = eval(XOUT.y);
    Yaw_pred = eval(XOUT.z);
    
    Updated_state.X_pred = X_pred;
    Updated_state.Y_pred = Y_pred;
    Updated_state.Yaw_pred = Yaw_pred;
end


