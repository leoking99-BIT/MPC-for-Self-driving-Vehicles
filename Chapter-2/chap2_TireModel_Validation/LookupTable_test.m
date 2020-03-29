function  LookupTable_test(u)

Alpha_L1    = u.signals.values(:,14);
Alpha_L2    = u.signals.values(:,15);

Fy_L1       = u.signals.values(:,22);
Fy_L2       = u.signals.values(:,23);
Num         = length(Fy_L2);

Fy_Alpha_ref = csvread('MFlookuptable_mu1_m1540.csv');

for i=1:1:Num
    Alpha_LT1(i) = func_LookupTable_alpha(Fy_Alpha_ref, Fy_L1(i)) ;
    Alpha_LT2(i) = func_LookupTable_alpha(Fy_Alpha_ref, Fy_L2(i)) ;
end

figure (1)
plot(1:Num, Alpha_LT1,'b',1:Num, Alpha_L1,'b*',1:Num, Alpha_L2,'r',1:Num, Alpha_LT2,'r*');
grid on

end


function [alpha] = func_LookupTable_alpha(Fy_Alpha_ref, Fy)
%***************************************************************%
% 1. Lookup a pre-defined table to convert the optimal Fyf to steering angle 
% 2. then transform the steering angle to Steer_SW
%***************************************************************%
	alpha = 0;
%     Fy_Alpha_ref = csvread('lookuptable_mu055_m1540.csv');
    TableSize = size(Fy_Alpha_ref); % row and column
    NumRow = TableSize(1,1);        %get the number of rows
    Index = 1;     
    while 1
        if Fy > Fy_Alpha_ref(1,1)
            alpha = Fy_Alpha_ref(1,2);
            break; 
        elseif Fy < Fy_Alpha_ref(NumRow,1)
            alpha = Fy_Alpha_ref(NumRow,2);
            break; 
        else % FyÔÚtable·¶Î§Ö®ÄÚ
            if Fy < Fy_Alpha_ref(Index,1)
                Index = Index+1;
            else
                k_slope = (Fy_Alpha_ref(Index,2)-Fy_Alpha_ref(Index-1,2))/(Fy_Alpha_ref(Index,1)-Fy_Alpha_ref(Index-1,1));
                alpha=Fy_Alpha_ref(Index-1,2)+(Fy-Fy_Alpha_ref(Index-1,1))* k_slope;              
                break;              
            end  
        end        
    end % end of while
%     alpha  = alpha*pi/180; %unit: deg-->rad
end