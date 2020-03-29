function [Updated_state] = func_UpdateState_EulerM_2_7(Previous_States, lfr, Vx_m_s, u, Ts)
% validation of dynamic model (2.26) use dsolve
    
    X_init      = Previous_States.X_pred;
    Y_init      = Previous_States.Y_pred;
    Yaw_init    = Previous_States.Yaw_pred;

    Updated_state.X_pred = X_init + Ts * Vx_m_s*cos(Yaw_init);
    Updated_state.Y_pred = Y_init + Ts * Vx_m_s*sin(Yaw_init);
    Updated_state.Yaw_pred = Yaw_init + Ts * Vx_m_s*tan(u)/lfr;
end