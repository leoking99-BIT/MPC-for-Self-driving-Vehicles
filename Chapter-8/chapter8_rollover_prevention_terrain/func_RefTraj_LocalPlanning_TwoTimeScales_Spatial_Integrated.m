function [WPIndex, RefP, RefU, Uaug, Uaug_0, PrjP, Roll_BaknR] = func_RefTraj_LocalPlanning_TwoTimeScales_Spatial_Integrated( MPCParameters, VehiclePara, WayPoints_Index, WayPoints_Collect, VehStateMeasured)
%***************************************************************%
%  不再有PrjP。
% 首先找到全局路径上距离车辆最近的点
% 其次，根据不同的步长，按照s选择一些列参考点并转换到车辆坐标系下。参考点的信息包括[s,x,y,bank]
% 再次，对车体坐标系下的x,y用Bezier曲线插值，定距采样，并计算采样点的航向角和曲率。4
% 同时，对道路的侧倾对s进行曲线拟合及采样
% 最后，将参考点的参数赋予RefP, RefU和Uaug

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
    Ns      = MPCParameters.Ns; % Tsplit
    Ts      = MPCParameters.Ts; % Set the sample time of near term
    Tsl     = MPCParameters.Tsl;% Set the sample time of long term   

    %------Measured or Estimated vehicle status
    Vel     = VehStateMeasured.x_dot;   % 20; % 
    PosX    = VehStateMeasured.X;
    PosY    = VehStateMeasured.Y;
    PosPsi  = VehStateMeasured.phi;      
    Roll_Shad = VehStateMeasured.Roll_Shad;%rad
%     Ax      = VehStateMeasured.Ax;
%     fwa     = VehStateMeasured.fwa; 
    
%*********** WaypointData2VehicleCoords ************************% 
    ds          = 0.1;%m
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
    %% 选择投影点， 全局路径上选择参考点，并转换到车体坐标下。同时选择s和对应的Bank angle        
        %--------------通过垂线找投影点--------------------------%
        [PPx,PPy,ey]=func_GetProjectPoint(WayPoints_Collect(index_min,2),... 
                                            WayPoints_Collect(index_min,3),... 
                                            WayPoints_Collect(index_min+1,2),... 
                                            WayPoints_Collect(index_min+1,3),... 
                                            PosX,... 
                                            PosY);
        Dy          = WayPoints_Collect(index_min+1,3) - WayPoints_Collect(index_min,3);
        Dx          = WayPoints_Collect(index_min+1,2) - WayPoints_Collect(index_min,2);
        Psi0        = atan2(Dy, Dx);  
        epsi        = Psi0 - PosPsi;
        
        PrjP.ey     = -ey;
        PrjP.epsi   = -epsi;        
        PrjP.Velr   = Vel;                                        
        PrjP.xr     = PPx;
        PrjP.yr     = PPy;
        PrjP.psir   = Psi0;
%         Kprj        = 0; 默认 投影点的曲率为零，即为一条直线
        PrjP.fwar   = 0; %atan(Kprj*L);  

        %-------------------i=1:Ns--根据车速在全局参考路径上选择参考点-------%
        Global_x        = [];
        Global_y        = [];  % 在全局路径上选择的路径点      
        Local_Sx        = [];
        Local_Sy        = [];
        Local_SS        = [];
        Local_SB        = []; %针对短时域
        StepLength_S    = Vel * Ts *  (Ns+1);% 多加一个点为了求曲线曲率时准备
        Ns_index        = index_min; 
        
%         %--index_min 与 index_min+1 之间的投影点
%         deltax          = PPx - PosX;
%         deltay          = PPy - PosY; 
%         CarCoord_x      = deltax * cos(PosPsi) + deltay * sin(PosPsi);
%         CarCoord_y      = deltay * cos(PosPsi) - deltax * sin(PosPsi);  
%         Local_Sx        = [Local_Sx; CarCoord_x];
%         Local_Sy        = [Local_Sy; CarCoord_y];           
%         Global_x        = [Global_x; PPx];
%         Global_y        = [Global_y; PPy];
%         
%         Local_SS        = [Local_SS; WayPoints_Collect(i,7)];
%         Local_SB        = [Local_SB; WayPoints_Collect(i,8)];       
%         %--index_min+1 
%         deltax          = WayPoints_Collect(index_min+1,2) - PosX;
%         deltay          = WayPoints_Collect(index_min+1,3) - PosY;
%         CarCoord_x      = deltax * cos(PosPsi) + deltay * sin(PosPsi);
%         CarCoord_y      = deltay * cos(PosPsi) - deltax * sin(PosPsi);    
%         Local_Sx        = [Local_Sx; CarCoord_x];
%         Local_Sy        = [Local_Sy; CarCoord_y];
%         Global_x        = [Global_x; WayPoints_Collect(index_min+1,2)];
%         Global_y        = [Global_y; WayPoints_Collect(index_min+1,3)];
%         
        tempDx          = WayPoints_Collect(index_min+1,2) - PPx;
        tempDy          = WayPoints_Collect(index_min+1,3) - PPy;
        Dist_1          = sqrt(power(tempDx,2) + power(tempDy,2)); %路点到投影点的距离 

        for i=index_min:1:WPNum %在参考路径上选择参考点,并通过坐标旋转转化到车体坐标系下
            Global_x        = [Global_x; WayPoints_Collect(i,2) ];
            Global_y        = [Global_y; WayPoints_Collect(i,3) ];  %先取出全局路径点
            
            deltax          = WayPoints_Collect(i,2) - PosX;
            deltay          = WayPoints_Collect(i,3) - PosY;
            CarCoord_x      = deltax * cos(PosPsi) + deltay * sin(PosPsi);
            CarCoord_y      = deltay * cos(PosPsi) - deltax * sin(PosPsi); % 全局路径点转换到局部坐标下              
            Local_Sx        = [Local_Sx; CarCoord_x];
            Local_Sy        = [Local_Sy; CarCoord_y];  %存储局部坐标下的点 
    
            
            Local_SS        = [Local_SS; WayPoints_Collect(i,7)];
            Local_SBL       = [Local_SB; WayPoints_Collect(i,8)];
            Local_SBR       = [Local_SB; WayPoints_Collect(i,9)];
            
            Ns_index        = i-1;            
            Dist_SumS       = Dist_1 + WayPoints_Collect(i,7) - WayPoints_Collect(index_min+1,7);  
            if(Dist_SumS >= StepLength_S)
                break;
            end            
        end % end of   for I=index_min+1:1:WPNum           
         
        %------------Ns:Np---------%
        Local_Lx        = [];
        Local_Ly        = [];
        Local_LS        = [];
        Local_LB        = [];
        StepLength_L    = Vel * Tsl * (Np-Ns+1);% 多加一个点为了求曲线曲率时准备
        Dist_SumL       = 0;      %
        for i=Ns_index:1:WPNum %在参考路径上选择参考点,并通过坐标旋转转化到车体坐标系下
            Global_x        = [Global_x; WayPoints_Collect(i,2) ];
            Global_y        = [Global_y; WayPoints_Collect(i,3) ];       
            
            deltax          = WayPoints_Collect(i,2) - PosX;
            deltay          = WayPoints_Collect(i,3) - PosY;
            CarCoord_x      = deltax * cos(PosPsi) + deltay * sin(PosPsi);
            CarCoord_y      = deltay * cos(PosPsi) - deltax * sin(PosPsi); % 转换到局部坐标下             
            Local_Lx        = [Local_Lx; CarCoord_x];
            Local_Ly        = [Local_Ly; CarCoord_y]; % 转换到局部坐标下 
  
            
            Local_LS        = [Local_LS; WayPoints_Collect(i,7)];
            Local_LBL       = [Local_LB; WayPoints_Collect(i,8)];
            Local_LBR       = [Local_LB; WayPoints_Collect(i,9)];
            
            Dist_SumL       = WayPoints_Collect(i,7) - WayPoints_Collect(Ns_index, 7 );
            if(Dist_SumL >= StepLength_L)
                break;
            end  
        end % end of   for i=Ns_index+1:1:WPNum   
        
        %%
        %------------多项式曲线拟合------------%
        if(Dist_SumS < StepLength_S) || (Dist_SumL < StepLength_L)
           WPIndex = 0; %如果没有找到则。。 % reaching the end ... %--这里没有考虑搜索到全局路径最后几个点时的情况，还不完备，有可能会报错！！！           
        else
             %----对短步长段Bezier曲线拟合，优点在于可以定距采样-----%
            MatS(:,1)=Local_Sx; 
            MatS(:,2)=Local_Sy;             
            [ps0,ps1,ps2,ps3,ts] = func_FindBezierControlPointsND(MatS,'u'); %uniform parameterization
            Scale                = round(Vel*Ts/ds);
            tlocS                = linspace(0,1,Scale*(Ns+1)+1);   %从起点到终点等距=0.1m采样,共（Ns+1）段，Scale*（Ns+1）+1个点
            MatLocalInterpS      = func_bezierInterp( ps0, ps1, ps2, ps3,tlocS);   % 曲线插值得到采样点
            
            MatSB(:,1)      = Local_SS; 
            MatSB(:,2)      = Local_SBR;             
            [psb0,psb1,psb2,psb3,tsb] = func_FindBezierControlPointsND(MatSB,'u'); %uniform parameterization
            tlocS                = linspace(0,1,Ns+2);   %从起点到终点等距采样,共（Np+1）段，（Np+2）个点
            MatLocalInterpSB     = func_bezierInterp( psb0,psb1,psb2,psb3,tlocS);   % 曲线插值得到采样点            
            
            Bezier_Sx       = zeros(Ns,1);
            Bezier_Sy       = zeros(Ns,1);
            Bezier_Spsi     = zeros(Ns,1);
            Bezier_SK       = zeros(Ns,1);
            Bezier_Sphi_t   = zeros(Ns,1);            
            for i = 2:1:length(MatLocalInterpSB(:,1))-1
                Bezier_Sx(i-1)    = MatLocalInterpS(Scale*(i-1),1);
                Bezier_Sy(i-1)    = MatLocalInterpS(Scale*(i-1),2);
                tempDx            = MatLocalInterpS(Scale*(i-1)+1,1) - MatLocalInterpS(Scale*(i-1),1);
                tempDy            = MatLocalInterpS(Scale*(i-1)+1,2) - MatLocalInterpS(Scale*(i-1),2);
                Bezier_Spsi(i-1)  = atan2(tempDy, tempDx);
                
%                 Bezier_SK(i-1)    = func_CalPathCurve(MatLocalInterpS(Scale*(i-1)-1,1),... 	% XA
%                                        MatLocalInterpS(Scale*(i-1)-1,2),...    % YA
%                                        MatLocalInterpS(Scale*(i-1),1),...      % XB
%                                        MatLocalInterpS(Scale*(i-1),2),...      % YB
%                                        MatLocalInterpS(Scale*(i-1)+1,1),...    % XC
%                                        MatLocalInterpS(Scale*(i-1)+1,2));      % YC
%  
                Bezier_SK(i-1)    = func_CalPathCurve_Patent(MatLocalInterpS(Scale*(i-1)-1,1),... 	% XA
                                       MatLocalInterpS(Scale*(i-1)-1,2),...    % YA
                                       MatLocalInterpS(Scale*(i-1),1),...      % XB
                                       MatLocalInterpS(Scale*(i-1),2),...      % YB
                                       MatLocalInterpS(Scale*(i-1)+1,1),...    % XC
                                       MatLocalInterpS(Scale*(i-1)+1,2));      % YC
%                 
%                 Bezier_SK(i-1)    = func_CalPathCurve_YU(MatLocalInterpS(Scale*(i-1)-1,1),... 	% XA
%                                        MatLocalInterpS(Scale*(i-1)-1,2),...    % YA
%                                        MatLocalInterpS(Scale*(i-1),1),...      % XB
%                                        MatLocalInterpS(Scale*(i-1),2),...      % YB
%                                        MatLocalInterpS(Scale*(i-1)+1,1),...    % XC
%                                        MatLocalInterpS(Scale*(i-1)+1,2));      % YC
                                   
                Bezier_Sphi_t(i-1) = MatLocalInterpSB(i,2);                  
            end % end of  for i = 2:1:length(MatLocalInterp(:,1))-1
            

            %----对长步长段 进行Bezier曲线拟合-----%
            MatL(:,1)=Local_Lx; 
            MatL(:,2)=Local_Ly;             
            [pL0,pL1,pL2,pL3,tL] = func_FindBezierControlPointsND(MatL,'u'); %uniform parameterization
            Scale                = round(Vel*Tsl/ds);
            tlocL                = linspace(0,1,Scale*(Np-Ns+1)+1);   %从起点到终点等距采样,共（Np-Ns+1）段，Scale*（Np-Ns+1）+1=Scale(Np-Ns)+11个点
            MatLocalInterpL      = func_bezierInterp( pL0, pL1, pL2, pL3,tlocL);   % 曲线插值得到采样点
            
            MatLB(:,1)=Local_LS; 
            MatLB(:,2)=Local_LBR;             
            [ps0,ps1,ps2,ps3,tL] = func_FindBezierControlPointsND(MatLB,'u'); %uniform parameterization
            tlocL                = linspace(0,1,Np-Ns+2);   %从起点到终点等距采样,共（Np+1）段，（Np+2）个点
            MatLocalInterpLB     = func_bezierInterp( ps0, ps1, ps2, ps3,tlocL);   % 曲线插值得到采样点    
            
            Bezier_Lx       = zeros(Np-Ns,1);
            Bezier_Ly       = zeros(Np-Ns,1);
            Bezier_Lpsi     = zeros(Np-Ns,1);
            Bezier_LK       = zeros(Np-Ns,1);
            Bezier_Lphi_t   = zeros(Np-Ns,1);  
            for i = 2:1:length(MatLocalInterpLB(:,1))-1
                Bezier_Lx(i-1)     = MatLocalInterpL(Scale*(i-1),1);
                Bezier_Ly(i-1)     = MatLocalInterpL(Scale*(i-1),2);
                tempDx             = MatLocalInterpL(Scale*(i-1)+1,1) - MatLocalInterpL(Scale*(i-1),1);
                tempDy             = MatLocalInterpL(Scale*(i-1)+1,2) - MatLocalInterpL(Scale*(i-1),2);
                Bezier_Lpsi(i-1)   = atan2(tempDy, tempDx);
                
%                 Bezier_LK(i-1)     = func_CalPathCurve(MatLocalInterpL(Scale*(i-1)-1,1),... 	% XA
%                                        MatLocalInterpL(Scale*(i-1)-1,2),...    % YA
%                                        MatLocalInterpL(Scale*(i-1),1),...      % XB
%                                        MatLocalInterpL(Scale*(i-1),2),...      % YB
%                                        MatLocalInterpL(Scale*(i-1)+1,1),...    % XC
%                                        MatLocalInterpL(Scale*(i-1)+1,2));      % YC
%                                    
                Bezier_LK(i-1)     = func_CalPathCurve_Patent(MatLocalInterpL(Scale*(i-1)-1,1),... 	% XA
                                       MatLocalInterpL(Scale*(i-1)-1,2),...    % YA
                                       MatLocalInterpL(Scale*(i-1),1),...      % XB
                                       MatLocalInterpL(Scale*(i-1),2),...      % YB
                                       MatLocalInterpL(Scale*(i-1)+1,1),...    % XC
                                       MatLocalInterpL(Scale*(i-1)+1,2));      % YC
%                                    
%                 Bezier_LK(i-1)     = func_CalPathCurve_YU(MatLocalInterpL(Scale*(i-1)-1,1),... 	% XA
%                                        MatLocalInterpL(Scale*(i-1)-1,2),...    % YA
%                                        MatLocalInterpL(Scale*(i-1),1),...      % XB
%                                        MatLocalInterpL(Scale*(i-1),2),...      % YB
%                                        MatLocalInterpL(Scale*(i-1)+1,1),...    % XC
%                                        MatLocalInterpL(Scale*(i-1)+1,2));      % YC                   
                                   
                Bezier_Lphi_t(i-1) = MatLocalInterpLB(i,2);
                
            end % end of  for i = 2:1:length(MatLocalInterp(:,1))-1


            
    %%
        RefP    = cell(Np,1);        
        RefU    = cell(Np,1); 
        for i = 1:1:Np
            if i <= Ns
               RefU{i,1} = atan(Bezier_SK(i)*L);                 
               RefP{i,1} = [Bezier_Sx(i);
                            Bezier_Sy(i);
                            Bezier_Spsi(i)]; 
            else
               RefU{i,1} = atan(Bezier_LK(i-Ns)*L); 
               RefP{i,1} = [Bezier_Lx(i-Ns);
                            Bezier_Ly(i-Ns);
                            Bezier_Lpsi(i-Ns)];                     
            end
        end  
        

        Uaug    = cell(Np,1);  
        
        Uaug_0  = [0;0]; % [MatLocalInterpSB(1,2); 0];  %          
        for i = 1:1:Np %不考虑道路曲率和倾角的影响
            if i <= Ns                     
                Uaug{i,1} = [0;0]; 
            else                   
                Uaug{i,1} = [0;0]; 
            end
        end  

%         Uaug_0  = [0; Bezier_SK(1)];  %  [0;0]; %      
%         for i = 1:1:Np %只考虑道路曲率
%             if i <= Ns                     
%                 Uaug{i,1} = [0; Bezier_SK(i)];
%             else                   
%                 Uaug{i,1} = [0; Bezier_LK(i-Ns)];
%             end
%         end  

%         Uaug_0  = [Roll_Shad; 0];  %  [MatLocalInterpSB(1,2); 0];  %  
%         for i = 1:1:Np %只考虑路面倾角的影响
%             if i <= Ns                     
%                 Uaug{i,1} = [Bezier_Sphi_t(i); 0];
%             else                   
%                 Uaug{i,1} = [Bezier_Lphi_t(i-Ns); 0];
%             end
%         end  
        
%         Uaug_0  = [Roll_Shad; Bezier_SK(1)];  %[MatLocalInterpSB(1,2); Bezier_SK(1)];  %      
%         for i = 1:1:Np %考虑道路曲率和倾角的影响
%             if i <= Ns                     
%                 Uaug{i,1} = [Bezier_Sphi_t(i); Bezier_SK(i)];
%             else                   
%                 Uaug{i,1} = [Bezier_Lphi_t(i-Ns); Bezier_LK(i-Ns)];
%             end
%         end  

        Roll_BaknR = MatLocalInterpSB(1,2);
        end % end of if(Dist_SumS < StepLength_S) || (Dist_SumL < StepLength_L)
        
    end % end of if( WPIndex > 0 )   % 如果找到了最近点

% %--------Plot local points and the fitted polynomial----------------%
% figure
% plot(Local_Sx,Local_Sy,'b*');
% hold on
% plot(Local_Lx,Local_Ly,'bo');    
% plot(Global_x,Global_y,'k.');  
% plot(Bezier_Sx,Bezier_Sy,'r+'); 
% plot(Bezier_Lx,Bezier_Ly,'ro'); 

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
% de与点到直线的距离不同，符号相反。
% 点到直线的距离：左正右负
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

function K=func_CalPathCurve(XA,YA,XB,YB,XC,YC)
    %% 通过三个点求圆心，再求圆心到任意一个点的距离
    %分别求两段直线的斜率
    if XB==XA
        mr=inf;
    else
        mr=(YB-YA)/(XB-XA);    
    end
    if XC==XB
        mt=inf;
    else
        mt=(YC-YB)/(XC-XB);    
    end

   %根据不同的斜率情况处理,mtdao=1/mt;mrsubmt=1/(2*(mr-mt));
    if mr==mt
        Rff=inf;
    else if mt==0
            if mr==inf
                Rff=sqrt(power((XA-XC),2)+power((YA-YC),2))/2;
            else
                mrsubmt=1/(2*(mr-mt));
                Xff=(mr*mt*(YC-YA)+mr*(XB+XC)-mt*mrsubmt*(XA+XB));
                Yff=(YB+YA)/2-(Xff-(XB+XA)/2)/mr;
                Rff=sqrt(power((XA-Xff),2)+power((YA-Yff),2));           
            end
        elseif mt==inf
            if mr==0
                Rff=sqrt(power((XA-XC),2)+power((YA-YC),2))/2;
            else
                Yff=(YB+YC)/2;
                Xff=(XA+XB)/2-mr*(YC-YA)/2;
                Rff=sqrt(power((XA-Xff),2)+power((YA-Yff),2));     
            end        
        else
            mtdao=1/mt;
            if mr==0
                mrsubmt=1/(2*(mr-mt));
                Xff=(mr*mt*(YC-YA)+mr*(XB+XC)-mt*mrsubmt*(XA+XB));
                Yff=(YB+YC)/2-mtdao*(Xff-(XB+XC)/2);
                Rff=sqrt(power((XA-Xff),2)+power((YA-Yff),2));            
            elseif mr==inf
                Yff=(YA+YB)/2;
                Xff=(XB+XC)/2+mt*(YA-YC)/2;
                Rff=sqrt(power((XA-Xff),2)+power((YA-Yff),2));
            else
                mrsubmt=1/(2*(mr-mt));
                Xff=(mr*mt*(YC-YA)+mr*(XB+XC)-mt*mrsubmt*(XA+XB));
                Yff=(YB+YC)/2-mtdao*(Xff-(XB+XC)/2);
                Rff=sqrt(power((XA-Xff),2)+power((YA-Yff),2));  
            end
        end    
    end
    %%
    %通过判断航向角的变化趋势来判定曲率的正负
    K1=GetPathHeading(XA,YA,XB,YB);
    K2=GetPathHeading(XB,YB,XC,YC);
    if K2>K1 %夹角变大，逆时针
        Rff=Rff;
    else %夹角变小，顺时针，K2<K1
        Rff=-1*Rff;
    end
    
    K = 1/Rff;
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

function curvature = func_CalPathCurve_YU(X1,Y1,X2,Y2,X3,Y3)
    %-----------caculate the radius of the circle first 
    % side one  
    delta_x = X2 - X1;  
    delta_y = Y2 - Y1;  
    a = sqrt(power(delta_x, 2.0) + power(delta_y, 2.0));  

    % side two  
    delta_x = X3 - X2;  
    delta_y = Y3 - Y2;  
    b = sqrt(power(delta_x, 2.0) + power(delta_y, 2.0));  

    % side three  
    delta_x = X1 - X3;  
    delta_y = Y1 - Y3;  
    c = sqrt(power(delta_x, 2.0) + power(delta_y, 2.0));  
    CLOSE_TO_ZERO = 0.01;
    if (a < CLOSE_TO_ZERO || b < CLOSE_TO_ZERO || c < CLOSE_TO_ZERO)   
        curvature = 0; 
    end
    
    %------------------------semiperimeter
    s = (a + b + c) / 2.0;  
    K = sqrt(abs(s * (s - a) * (s - b) * (s - c)));  
    curvature = 4 * K / (a * b * c);  
    
    %------------ determine the sign, using cross product(叉乘)
    % 2维空间中的叉乘是： A x B = |A||B|Sin(θ)
    % V1(x1, y1) X V2(x2, y2) = x1y2 – y1x2  
    rotate_direction = (X2 - X1) *  (Y3 - Y2) - (Y2 - Y1) * (X3 - X2);  
    if(rotate_direction < 0) %通过判断符号，判定Sin(θ)的符号
        %Sin(θ)<0, 顺时针旋转，曲率为负
        curvature = -curvature;  
    end

end








