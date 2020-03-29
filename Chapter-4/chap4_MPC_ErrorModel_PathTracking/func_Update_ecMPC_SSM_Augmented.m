function [StateSpaceModel] = func_Update_ecMPC_SSM_Augmented(VehiclePara, MPCParameters, Vel)
    % generate State-space model
   

lfr = VehiclePara.L; 
Ts  = MPCParameters.Ts;

Au = [1,        0,  Vel*Ts/lfr;
      -Vel*Ts,  1,  0;
      0,        0,  1];
Bu1 = [Vel*Ts/lfr; 0; 1];
Bu2 = [-Vel*Ts/lfr; 0; 0];
C = [1, 0, 0;
     0, 1, 0];
 
    StateSpaceModel.Au  = Au;
    StateSpaceModel.Bu1 = Bu1; 
    StateSpaceModel.Bu2 = Bu2;
    StateSpaceModel.C   = C;    
end % EoF

