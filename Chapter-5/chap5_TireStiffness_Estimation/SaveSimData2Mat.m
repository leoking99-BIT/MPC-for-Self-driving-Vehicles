
% Sim_Data_Both_ConstantCorneringStiff = Sim_Data_Both_ConstantCornering;

% save SimResult_Both.mat Sim_Data;

% save SimResult_OnlyBank.mat Sim_Data;

% save SimResult_OnlyCurvature.mat Sim_Data;
% 

Outputed_Data = Outputed_Data.signals.values(:,:);

save TireStiffnessEstimation_Data.mat Outputed_Data;

