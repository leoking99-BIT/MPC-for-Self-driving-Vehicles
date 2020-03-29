function [Envelope] = func_SafedrivingEnvelope_SL(VehiclePara, MPCParameters, Constraints, StateSpaceModel, Vel, CarHat, EY_MAX, EY_Min)
    %  Generating safe driving envelope 
    h       = VehiclePara.hCG;
    Tr      = VehiclePara.Tr;
    lf      = VehiclePara.Lf;
    lr      = VehiclePara.Lr;
    l       = VehiclePara.L;
    M       = VehiclePara.m;       
    g       = VehiclePara.g;
    mu      = VehiclePara.mu;
    Iz      = VehiclePara.Iz;   %I为车辆绕Z轴的转动惯量，车辆固有参数  
    Ix      = VehiclePara.Ix;   %I为车辆绕Z轴的转动惯量，车辆固有参数  
    Np      = MPCParameters.Np;
    Ns      = MPCParameters.Ns;
    Alphar_lim = Constraints.arlim;
    
    %-----------environmental envelope-------------%
    Envelope.Henv  = [0   0   0   0   1     0;
                      0   0   0   0   -1    0]; 
                  
    Envelope.Genv  = cell(Np,1); 
    Dist_temp      = 0.8; % dw/2 + dbuffer; 
    for i = 1:1:Np
       Envelope.Genv{i,1} = [EY_MAX(i) - Dist_temp;
                             -EY_Min(i) - Dist_temp]; 
    end  

    %-----------stable handling envelope-------------%
%     for j = 1:1:Np
        Envelope.Hsh  = [1/Vel   -lr/Vel    0   0    0   0;
                         0         1        0   0    0   0];   
        Envelope.Psh  = [0  0;  g/Vel   0];
        
        r_ssmax = -CarHat*Alphar_lim*(1+lr/lf)/(M*Vel);                  

        Envelope.Gsh  = [Alphar_lim    r_ssmax]';                  
%     end

   %-----------zero moment point -------------% 
   Acn = StateSpaceModel.Acn;  % 8*8
   B1cn = StateSpaceModel.B1cn; %8*1
   B2cn = StateSpaceModel.B2cn;%8*2
   
%    N1 = [h/g    0            -Ix/(M*g)   0   0   0];
%    N2 = [0      h*Vel/g      0           h   0   0]; %
%    N3 = [h  0];

   N1 = [2*h/(g*Tr)     0                   -2*(M*h*h+Ix)/(M*g*Tr)      0         0   0];
   N2 = [0              2*h*Vel/(g*Tr)     0                            2*h/Tr    0   0]; %
%    N3 = [2*h/Tr         0];

   
   H_yzmp   = N1*Acn + N2;%1*8
   P_yzmp1  = N1*B1cn; %1*1
   P_yzmp2  = N1*B2cn; % + N3;%1*2
   
   Envelope.H_yzmp  = H_yzmp;
   Envelope.P_yzmp1 = P_yzmp1;
   Envelope.P_yzmp2 = P_yzmp2;
    
end  % end of func_SafedrivingEnvelope 