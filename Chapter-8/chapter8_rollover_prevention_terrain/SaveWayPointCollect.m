
% WayPoints_Collect = u.signals.values(1:2600,:);
% vel, x,   y,  psi, ax,   steering, station, BankL, BankR
% save WayPoints_Alt3fromFHWA_Overall_Station_Bank.mat WayPoints_Collect;

% MatL = MatL(1:21,:);
% save Bankangle_LeftLane_Alt3fromFHWA.mat MatL;
% 
% MatR = MatR(1:21,:);
% save eBankangle_RightLan_Alt3fromFHWA.mat MatR;


% DataNoConstr =u.signals.values(2:601,:);
% save CurveTracking_NoConstr_20.mat  DataNoConstr

DataWithConstr =u.signals.values(2:601,:);
save CurveTracking_WithConstr_20.mat  DataWithConstr