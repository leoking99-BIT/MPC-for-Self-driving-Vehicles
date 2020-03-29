function [ U1_Index, U1] = func_LocalControlInput_TwoTimeScales(MPCParameters, VehiclePara, index, fwa)

%*********** Parameters Initialization *************************% 
    L       = VehiclePara.L;   % ≥µ¡æ÷·æ‡
    Np      = MPCParameters.Np;% ‘§≤‚ ±”Ú
    Ns      = MPCParameters.Ns; % Tsplit
    Ts      = MPCParameters.Ts; % Set the sample time of near term
    Tsl     = MPCParameters.Tsl;% Set the sample time of long term   

%     U1_0    = fwa(index+1);
%     TScale  = Tsl/Ts;
%     U1      = zeros(Np+1,1);
%     U1_Index = zeros(Np+1,1);
%     for i = 1:1:Np+1
%         if i <= Ns  % u0~u_{Ns-1}
%             U1_Index(i)  = index+i;
%             U1(i)       = fwa(index+i);
%         else
%             temp        = TScale*(i-Ns)+Ns;
%             U1_Index(i)  = index+temp;
%             U1(i)       = fwa(index+temp);     
%         end
%     end  

    U1          = zeros(Np+1,1);
    U1_Index    = zeros(Np+1,1);
    for i = 1:1:Np+1
        U1_Index(i)  = index+i;
        U1(i)        = fwa(index+i);        
    end

    
end