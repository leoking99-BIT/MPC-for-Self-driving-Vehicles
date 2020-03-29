function Generate1DLookuptable_BT()
%------------------------------------------------------------------%
% Brush Tire Model by Pacejka
% Using Brush tire model to calcualte Fy 
% Assume slip = 0, 
%------------------------------------------------------------------%
    mu = 1.0;  %地面摩擦因数
%********(1)calculate normal tire load **************%
    m  = 1540; %为车辆质量,Kg
    g  = 9.8;
    lf = 1.11; %unit:m
    lr = 1.67; %unit:m 前后车轮距离车辆质心的距离，车辆固有参数
    L = 2.78;  %unit:m
    Fzf = 0.5*m*g*lr/L;   
%********(2)设定参数 **************%
    Fz0=4100;                     %nominal (rated) load(>0）,N  %额定垂直载荷
    Pky1=-12.95;Pky2=1.72;Pky3=0.22; 
    r=0; %0.1*pi/180;                 %camber angle 侧倾角
    sr=sin(r);                        %r 是外倾角，sr 表示r*
    Kya0=Pky1*Fz0*sin(2*atan(Fzf/(Pky2*Fz0)));%取lam（Kya）=1 
    Calpha = -1*Kya0*(1-Pky3*sr^2);  
    Ca = 2* Calpha;
%********(3)alpha processing **************%
%     Alpha_max = atan(3*1540*9.8*0.55*1.11/(140000*2.78))*180/pi; % 4deg
%     alphamax = atan(3*mu*Fzf/Ca)*180/pi; % 约等于8deg
%     alpha = -alphamax:0.1:alphamax;
    alpha = -8:0.1:8;
    Alpha_rad = alpha.*pi/180;% % alpha  in deg, transform into rad
    tanalpha_rad = tan(Alpha_rad);
%********(4) Using Brush tire model to calcualte Fy  **************%
    Fyf = -Calpha * tanalpha_rad + power(Calpha,2) * tanalpha_rad.*abs(tanalpha_rad)/(3*mu*Fzf) - power(Calpha,3)*power(tanalpha_rad,3)/(27*mu*mu*Fzf*Fzf); 
    FyBT  = Fyf;

%     table(:,1) = FyBT;
%     table(:,2) = alpha;
%     csvwrite('BTlookuptable_mu1_m1540.csv',table); 

%********(5) Using MF tire model to calcualte Fy  **************% 
num = length(alpha);
FyMF = zeros(num,1);
for i = 1:1:num
    FyMF(i) = func_MFTyreModel_puresideslip(Fzf, alpha(i));
end

    plot(alpha, FyBT,'b',alpha, FyMF,'r');
    grid on

%     table(:,1) = FyMF;
%     table(:,2) = alpha;
%     csvwrite('MFlookuptable_mu1_m1540.csv',table);     
    
%     Alpha_max = atan(3*1973*9.8*0.55*1.53/(140000*2.76))*180/pi;
%     Alpha_max = atan(3*1540*9.8*0.55*1.11/(140000*2.78))*180/pi;

    
end