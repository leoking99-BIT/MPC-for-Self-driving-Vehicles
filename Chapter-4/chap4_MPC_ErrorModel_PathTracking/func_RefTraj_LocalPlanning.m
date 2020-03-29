function [WPIndex, RefP, RefK, Uaug, PrjP] = func_RefTraj_LocalPlanning( MPCParameters, VehiclePara, WayPoints_Index, WayPoints_Collect, VehStateMeasured)
%***************************************************************%
% 首先找到全局路径上距离车辆最近的点 (相当于投影点)
% 其次，根据步长，按照s选择一些列参考点并转换到车辆坐标系下。参考点的信息包括[s,x,y]
% 再次，对车体坐标系下的x,y用Bezier曲线插值，定距采样，并计算采样点的航向角和曲率。
% 最后，将参考点的参数赋予RefP, RefU和Uaug
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
% Email:leoking1025@gmail.com
% My homepage: https://sites.google.com/site/kailiumiracle/  
%***************************************************************%

%*********** Parameters Initialization *************************% 
    L       = VehiclePara.L;   % 车辆轴距
    Np      = MPCParameters.Np;% 预测时域
    Ts      = MPCParameters.Ts; % Set the sample time

    %------Measured or Estimated vehicle status
    Vel     = VehStateMeasured.x_dot;
    PosX    = VehStateMeasured.X;
    PosY    = VehStateMeasured.Y;
    PosPsi  = VehStateMeasured.phi;      
    
%*********** WaypointData2VehicleCoords ************************% 
    ds          = 0.1;%unit:m, 路径点之间的距离
    WPNum       = length(WayPoints_Collect(:,1));
    
    %--------先找到参考路径上距离车辆最近的点--------------------------%  
    Dist_MIN    = 10000;
    index_min   = 0;
    for i = WayPoints_Index:1:WPNum 
        deltax  = WayPoints_Collect(i,2) - PosX;
        deltay  = WayPoints_Collect(i,3) - PosY;
        Dist    = sqrt(power(deltax,2) + power(deltay,2));% 路点到车辆重心的距离
        if Dist < Dist_MIN
            Dist_MIN  = Dist; 
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
    %% 选择投影点， 全局路径上选择参考点，并转换到车体坐标下。   
        %--------------通过垂线找投影点--------------------------%
        [PPx,PPy,ey]=func_GetProjectPoint(WayPoints_Collect(index_min,2),... 
                                            WayPoints_Collect(index_min,3),... 
                                            WayPoints_Collect(index_min+1,2),... 
                                            WayPoints_Collect(index_min+1,3),... 
                                            PosX,... 
                                            PosY);
        Dy          = WayPoints_Collect(index_min+1,3) - WayPoints_Collect(index_min,3);
        Dx          = WayPoints_Collect(index_min+1,2) - WayPoints_Collect(index_min,2);
        Psi0        = atan2(Dy, Dx); % [-pi, pi]
        epsi        = PosPsi - Psi0;% 车辆航向 - 道路切向，（逆时针为正）
        
        PrjP.epsi   = epsi;%
        PrjP.ey     = ey;%ey的方向定义为左负右正,车辆相对于参考道路
        PrjP.Velr   = Vel;                                        
        PrjP.xr     = PPx;
        PrjP.yr     = PPy;
        PrjP.psir   = Psi0;
        PrjP.fwar   = 0; %atan(Kprj*L);  

        %-------------------i=1:Np--根据车速在全局参考路径上选择参考点-------%
        Local_Sx        = [];
        Local_Sy        = [];
        StepLength_S    = Vel * Ts *  (Np+1);% 多加一个点为了求曲线曲率时准备
            
        tempDx          = WayPoints_Collect(index_min+1,2) - PPx;
        tempDy          = WayPoints_Collect(index_min+1,3) - PPy;
        Dist_1          = sqrt(power(tempDx,2) + power(tempDy,2)); %路点到投影点的距离 

        for i=index_min:1:WPNum %在参考路径上选择参考点,并通过坐标旋转转化到车体坐标系下
            deltax          = WayPoints_Collect(i,2) - PosX;
            deltay          = WayPoints_Collect(i,3) - PosY;
            CarCoord_x      = deltax * cos(PosPsi) + deltay * sin(PosPsi);
            CarCoord_y      = deltay * cos(PosPsi) - deltax * sin(PosPsi); % 全局路径点转换到局部坐标下              
            Local_Sx        = [Local_Sx; CarCoord_x];
            Local_Sy        = [Local_Sy; CarCoord_y];  %存储局部坐标下的点  
                    
            Dist_SumS       = Dist_1 + WayPoints_Collect(i,7) - WayPoints_Collect(index_min+1,7);  
            if(Dist_SumS >= StepLength_S)
                break;
            end            
        end % end of   for I=index_min+1:1:WPNum           
        
        %%
        %------------多项式曲线拟合------------%
        if(Dist_SumS < StepLength_S)
           WPIndex = 0; %如果没有找到则。。 % reaching the end ... %--这里没有考虑搜索到全局路径最后几个点时的情况，还不完备，有可能会报错！！！           
        else
             %----Bezier曲线拟合，优点在于可以定距采样-----%
            MatS(:,1)=Local_Sx; 
            MatS(:,2)=Local_Sy;             
            [ps0,ps1,ps2,ps3,ts] = func_FindBezierControlPointsND(MatS,'u'); %uniform parameterization
            Scale                = round(Vel*Ts/ds);
            tlocS                = linspace(0,1,Scale*(Np+1)+1);   %从起点到终点等距=0.1m采样,共（Np+1）段，Scale*（Np+1）+1个点
            MatLocalInterpS      = func_bezierInterp( ps0, ps1, ps2, ps3,tlocS);   % 曲线插值得到采样点       
            
            Bezier_Sx       = zeros(Np,1);
            Bezier_Sy       = zeros(Np,1);
            Bezier_Spsi     = zeros(Np,1);
            Bezier_SK       = zeros(Np,1);         
            for i = 2:1:Np+1
                Bezier_Sx(i-1)    = MatLocalInterpS(Scale*(i-1),1);
                Bezier_Sy(i-1)    = MatLocalInterpS(Scale*(i-1),2);
                tempDx            = MatLocalInterpS(Scale*(i-1)+1,1) - MatLocalInterpS(Scale*(i-1),1);
                tempDy            = MatLocalInterpS(Scale*(i-1)+1,2) - MatLocalInterpS(Scale*(i-1),2);
                Bezier_Spsi(i-1)  = atan2(tempDy, tempDx);
                
                Bezier_SK(i-1)    = func_CalPathCurve_Patent(MatLocalInterpS(Scale*(i-1)-1,1),... 	% XA
                                       MatLocalInterpS(Scale*(i-1)-1,2),...    % YA
                                       MatLocalInterpS(Scale*(i-1),1),...      % XB
                                       MatLocalInterpS(Scale*(i-1),2),...      % YB
                                       MatLocalInterpS(Scale*(i-1)+1,1),...    % XC
                                       MatLocalInterpS(Scale*(i-1)+1,2));      % YC            
            end % end of  for i = 2:1:length(MatLocalInterp(:,1))-1
            
    %%
        RefP    = cell(Np,1);        
        RefK    = cell(Np,1); 
        Uaug    = cell(Np,1);      
        for i = 1:1:Np
           Uaug{i,1} = atan(Bezier_SK(i)*L);  
           RefK{i,1} = -Bezier_SK(i);     
           RefP{i,1} = [Bezier_Sx(i);
                        Bezier_Sy(i);
                        Bezier_Spsi(i)]; 
        end

        end % end of if(Dist_SumS < StepLength_S) || (Dist_SumL < StepLength_L)
        
    end % end of if( WPIndex > 0 )   % 如果找到了最近点

end % end of function 


%==============================================================%
% sub functions
%==============================================================%   
function K=GetPathHeading(Xb,Yb,Xn,Yn)
    %***Way I.求Heading Angle 在[-pi,pi]之间 *******%
    AngleY=Yn-Yb;
    AngleX=Xn-Xb;
    K= atan2(AngleY, AngleX);
    
    %***Way II. 求Heading Angle 在0~2*pi之间 *******%
%     AngleY=Yn-Yb;
%     AngleX=Xn-Xb;    
%     
%     if Xb==Xn
%         if Yn>Yb
%             K=pi/2;
%         else
%             K=3*pi/2;
%         end
%     else
%         if Yb==Yn
%             if Xn>Xb
%                 K=0;
%             else
%                 K=pi;
%             end
%         else
%             K=atan(AngleY/AngleX);
%         end    
%     end
% 
%     %****修正K,使之在0~360°之间*****%
%    if (AngleY>0&&AngleX>0)%第一象限
%         K=K;
%     elseif (AngleY>0&&AngleX<0)||(AngleY<0&&AngleX<0)%第二、三象限
%         K=K+pi;
%     else if (AngleY<0&&AngleX>0)%第四象限
%             K=K+2*pi;  
%         else
%             K=K;
%         end
%    end
    
end % end of function

function [PPx,PPy,de]=func_GetProjectPoint(Xb,Yb,Xn,Yn,Xc,Yc)
%-------------------------------------------------------%
% 通常点到直线的距离定义为：左正右负
% 这里de与点到直线的距离的方向定义不同，符号相反。
% 即， de的方向定义为左负右正
%-------------------------------------------------------%

    if Xn==Xb
        x=Xn;
        y=Yc;
        de=Xc-Xn;
    else if Yb==Yn
            x=Xc;
            y=Yn;
            de=Yn-Yc;
        else
            DifX=Xn-Xb;
            DifY=Yn-Yb;
            Kindex=DifY/DifX;
            bindex=Yn-Kindex*Xn;
            
            K=(-1)*1/Kindex;
            b=Yc-K*Xc;
            x=(bindex-b)/(K-Kindex);
            y=K*x+b;
            de=(Kindex*Xc+bindex-Yc)/sqrt(1+Kindex*Kindex);
        end     
    end
    PPx=x;
    PPy=y;
       
end

function K=func_CalPathCurve_Patent(XA,YA,XB,YB,XC,YC)
    %% 
    x_dot       = XC - XA;
    y_dot       = YC - YA;
    x_dotdot    = XC + XA - 2*XB;
    y_dotdot    = YC + YA - 2*YB;
    temp        = x_dot*x_dot + y_dot*y_dot;
    K= 4*(x_dot*y_dotdot - x_dotdot*y_dot )/ power(temp, 1.5);

end








