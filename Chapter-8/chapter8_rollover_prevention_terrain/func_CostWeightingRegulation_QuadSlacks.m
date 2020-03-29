function [Sl, Ql, Rdun, Wshl, dun, dul] = func_CostWeightingRegulation_QuadSlacks(MPCParameters, CostWeights, Constraints)

%% 参数初始化
    Ts      = MPCParameters.Ts;
    Tsl     = MPCParameters.Tsl;
    
    Qephi   = CostWeights.Wephi;
    Qey     = CostWeights.Wey;
    RDdeltaf= CostWeights.Ddeltaf;   
    Sdeltaf = CostWeights.deltaf;    
    Wshar   = CostWeights.Wshar;
    Wshr    = CostWeights.Wshr;

    DPhimax = Constraints.DPhimax;  %  0.15 rad ==> 8.5deg
    Dymax   = Constraints.Dymax;
    armax   = Constraints.arlim; % 0.104rad=6deg; %0.15rad=8deg
    rmax    = Constraints.rmax;    

    umax    = Constraints.umax;
    dumax   = Constraints.dumax;
    dun     = dumax * Ts;
    dul     = dumax * Tsl;
    
    %% 权重因子归一化
    Qephi_DPhimax2  = Qephi/(DPhimax*DPhimax);
    Qey_Dymax2      = Qey/(Dymax*Dymax);    
    Ql              = diag([0, 0, 0, 0, Qephi_DPhimax2, Qey_Dymax2]);
    
    Wshr_rmax2      = Wshr/(rmax*rmax);    
    Wshar_armax2    = Wshar/(armax*armax);
    Wshl            = diag([Wshar_armax2, Wshr_rmax2]);
   
    dumax_ts2   = (dumax * dumax* Ts * Ts);% 
    Rdun          = RDdeltaf/dumax_ts2;
%     Rdul        = Ts_Tsl * Rdun; 
    
    Sl = Sdeltaf/(umax * umax);
    
    
end % end of func_CostWeightingRegulation