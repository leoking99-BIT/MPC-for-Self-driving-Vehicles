
% Sim_Data_Both_ConstantCorneringStiff = Sim_Data_Both_ConstantCornering;

% save SimResult_Both.mat Sim_Data;

% save SimResult_OnlyBank.mat Sim_Data;

% save SimResult_OnlyCurvature.mat Sim_Data;
% 

mu_dot3_50 = u.signals.values(:,:);

save SimData_mu_dot3_50.mat mu_dot3_50;

