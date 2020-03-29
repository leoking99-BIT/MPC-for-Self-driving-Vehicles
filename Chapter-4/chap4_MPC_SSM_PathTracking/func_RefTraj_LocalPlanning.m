function [PrjP, RefP, RefU, WPIndex ] = func_RefTraj_LocalPlanning( MPCParameters, WayPoints_Index, WayPoints_Collect, VehStateMeasured)
%***************************************************************%
%  
% Input:
% MPCParameters：
% WayPoints_Index
% WayPoints_Collect
% VehStateMeasured
% 
% Output:
% WPIndex： 
%     > 0：Normal, WPIndex = index_min;
%     0:error,index_min<1
%     -1: index_min = WPNum,到了全局路径的尽头，停车
%---------------------------------------------------------------%
% Published by: Kai Liu
% Email:leoking1025@bit.edu.cn
% My github: https://github.com/leoking99-BIT 
%***************************************************************%

%*********** Parameters Initialization *************************% 
    L       = 2.78;
    Np      = MPCParameters.Np;
    Ts      = MPCParameters.Ts;
    %------Measured or Estimated vehicle status
    Vel     = VehStateMeasured.Vx;   % 20; % 
    PosX    = VehStateMeasured.Xc;
    PosY    = VehStateMeasured.Yc;
    PosPsi  = VehStateMeasured.psi;            
    Ax      = VehStateMeasured.Ax;
    fwa     = VehStateMeasured.fwa; 
    
    RefP    = cell(Np,1); 
    RefU    = cell(Np,1);   
%*********** WaypointData2VehicleCoords ************************% 
    WPNum       = length(WayPoints_Collect(:,1));
    
    %--------先找到参考路径上距离车辆最近的点--------------------------%  
    Dist_MIN    = 1000;
    index_min   = 0;
    for i=WayPoints_Index:1:WPNum 
        deltax  = WayPoints_Collect(i,2) - PosX;
        deltay  = WayPoints_Collect(i,3) - PosY;
        Dist    = sqrt(power(deltax,2) + power(deltay,2)); %路点到车辆重心的距离
        if Dist < Dist_MIN
            Dist_MIN = Dist; 
            index_min = i;
        end
    end
    if (index_min < 1) 
        WPIndex = -1; %如果没有找到则。。
    else if ( index_min >= WPNum)
            WPIndex = -2; %如果没有找到则。。
        else
            WPIndex = index_min;
        end
    end
    
    
    if( WPIndex > 0 )   % 如果找到了最近点
        %--------------通过垂线找投影点--------------------------%
        [PPx,PPy,ey]=func_GetProjectPoint(WayPoints_Collect(index_min,2),... 
                                            WayPoints_Collect(index_min,3),... 
                                            WayPoints_Collect(index_min+1,2),... 
                                            WayPoints_Collect(index_min+1,3),... 
                                            PosX,... 
                                            PosY);
        Psi0 = atan2(WayPoints_Collect(index_min+1,3) - PPy, WayPoints_Collect(index_min+1,2) - PPx);
        epsi        = Psi0 - PosPsi;
        PrjP.ey     = ey;
        PrjP.epsi   = epsi;
        PrjP.xr     = PPx;
        PrjP.yr     = PPy;
        PrjP.psir   = Psi0;
        PrjP.Velr   = Vel; 
        PrjP.fwar   = 0; %atan(Kprj*L);
        
        %------------根据车速在全局参考路径上选择参考点--------------------%
        StepLength      = Vel * Ts * (Np+1);% 多加一个点为了求曲线曲率时准备
        
        count           = 1; % 将投影点作为曲线的起点
        deltax          = PPx - PosX;
        deltay          = PPy - PosY; 
        Local_x(count)  = deltax * cos(PosPsi) + deltay * sin(PosPsi);
        Local_y(count)  = deltay * cos(PosPsi) - deltax * sin(PosPsi);   
        Dist_Sum        = 0;      %
        
        count           = count + 1; % 将投影点之后的第一个点转换到局部坐标下   
        deltax          = WayPoints_Collect(index_min+1,2) - PosX;
        deltay          = WayPoints_Collect(index_min+1,3) - PosY;
        Local_x(count)  = deltax * cos(PosPsi) + deltay * sin(PosPsi);
        Local_y(count)  = deltay * cos(PosPsi) - deltax * sin(PosPsi);          
        tempDx          = WayPoints_Collect(index_min+1,2) - PPx;
        tempDy          = WayPoints_Collect(index_min+1,3) - PPy;
        DistBtn         = sqrt(power(tempDx,2) + power(tempDy,2)); %路点到投影点的距离
        Dist_Sum        = Dist_Sum + DistBtn;            
        
        for i=index_min+2:1:WPNum %在参考路径上选择参考点,并通过坐标旋转转化到车体坐标系下
            count           = count + 1;            
            deltax          = WayPoints_Collect(i,2) - PosX;
            deltay          = WayPoints_Collect(i,3) - PosY;
            Local_x(count)  = deltax * cos(PosPsi) + deltay * sin(PosPsi);
            Local_y(count)  = deltay * cos(PosPsi) - deltax * sin(PosPsi); % 转换到局部坐标下   

            tempDx          = WayPoints_Collect(i,2) - WayPoints_Collect(i-1,2);
            tempDy          = WayPoints_Collect(i,3) - WayPoints_Collect(i-1,3);
            DistBtn         = sqrt(power(tempDx,2) + power(tempDy,2)); %此路点到前一个路点的距离
            Dist_Sum        = Dist_Sum + DistBtn;    
            
            if(Dist_Sum > StepLength)
                break;
            end            
        end % end of   for j=index_min+1:1:WPNum       
        
        %------------多项式曲线拟合------------%
        if(Dist_Sum < StepLength)
           WPIndex = 0; %如果没有找到则。。 % reaching the end ... %--这里没有考虑搜索到全局路径最后几个点时的情况，还不完备，有可能会报错！！！
        else %在全局路径上找到了合适的参考点，多项式拟合
            %----Bezier曲线拟合，优点在于可以定距采样-----%
            Mat(:,1)=Local_x'; 
            Mat(:,2)=Local_y'; 
            
            [p0,p1,p2,p3,t] = func_FindBezierControlPointsND(Mat,'u'); %uniform parameterization
            tloc            = linspace(0,1,Np+2);   %从起点到终点等距采样,共（Np+1）段，（Np+2）个点
            MatLocalInterp  = func_bezierInterp( p0, p1, p2, p3,tloc);   % 曲线插值得到采样点
            
            Bezier_x    = zeros(Np,1);
            Bezier_y    = zeros(Np,1);
            Bezier_psi  = zeros(Np,1);
            Bezier_K    = zeros(Np,1);
            for i = 2:1:length(MatLocalInterp(:,1))-1
                Bezier_x(i-1)     = MatLocalInterp(i,1);
                Bezier_y(i-1)     = MatLocalInterp(i,2);
                tempDx            = MatLocalInterp(i+1,1) - MatLocalInterp(i,1);
                tempDy            = MatLocalInterp(i+1,2) - MatLocalInterp(i,2);
                Bezier_psi(i-1)   = atan2(tempDy, tempDx);
                
                Bezier_K(i-1)     = func_CalPathCurve(MatLocalInterp(i-1,1),... 	% XA
                                       MatLocalInterp(i-1,2),...    % YA
                                       MatLocalInterp(i,1),...      % XB
                                       MatLocalInterp(i,2),...      % YB
                                       MatLocalInterp(i+1,1),...    % XC
                                       MatLocalInterp(i+1,2));      % YC
            end % end of  for i = 2:1:length(MatLocalInterp(:,1))-1
        end    % end of if(Dist_Sum < StepLength)       
 
        for i = 1:1:Np
           RefP{i,1} = [Bezier_x(i);
                        Bezier_y(i);
                        Bezier_psi(i)]; 
        end  

        for i = 1:1:Np
           RefU{i,1} = atan(Bezier_K(i)*L); 
        end  
    end % end of if (index_min < 1) || ( index_min > WPNum)


end % end of function 


%==============================================================%
% sub functions
%==============================================================%   

%--------Plot local points and the fitted polynomial----------------%
% Fitted_y = K(1)*Local_x.*Local_x.*Local_x + K(2)*Local_x.*Local_x + K(3)*Local_x + K(4);
% % y = K(4)*x.*x.*x+K(3)*x.*x+K(2)*x+K(1);
% figure(1)
% plot(Local_x,Local_y,'b')
% hold on
% plot(Local_x,Fitted_y,'r');





