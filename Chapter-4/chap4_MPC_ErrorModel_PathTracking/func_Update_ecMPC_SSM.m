%***************************************************************%
% Deem Vx as constant
%***************************************************************%
function [StateSpaceModel] = func_Update_ecMPC_SSM(VehiclePara, MPCParameters, Vel)
    % generate State-space model
   

lfr = VehiclePara.L;  %a = 1.11;  
Ts  = MPCParameters.Ts;

% Au = [1,        0,  Vel*Ts/lfr;
%       -Vel*Ts,  1,  0;
%       0,        0,  1];
% Bu1 = [Vel*Ts/lfr; 0; 1];
% Bu2 = [-Vel*Ts/lfr; 0; 0];
% C = [1, 0, 0;
%      0, 1, 0];

Au = [1,        0;
      -Vel*Ts,  1];
Bu1 = [Vel*Ts/lfr; 0];
Bu2 = [-Vel*Ts/lfr; 0];
 
    StateSpaceModel.Au  = Au;
    StateSpaceModel.Bu1 = Bu1; 
    StateSpaceModel.Bu2 = Bu2;
    
end % EoF

