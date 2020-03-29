function [Fy0] = func_MFTyreModel_puresideslip(Fz, alpha)
%------------------------------------------------------------------%
% Magic Formula model 
%------------------------------------------------------------------%
% Parameter values of Magic Formula from CarSim
Fz0=4100;                     %nominal (rated) load(>0）,N  %额定垂直载荷
R0=0.298;                     %unloaded tyre radius (=ro)，m  
dfz=(Fz-Fz0)/Fz0;  
mux = 1;  muy = 1;            % Tire/ground coefficients for this data
% Slip_Ratio = 0;             % pure side slip
%-------Pacejka Tire 175/70 R13(symmetric)----------------------%
Pcy1=1.29;  
Pdy1=-0.9;Pdy2=0.18;Pdy3=-4.5;  
Pey1=-1.07;Pey2=0.68;Pey3=-0.63;Pey4=-12.35;
Pky1=-12.95;Pky2=1.72;Pky3=0.22;     
Phy1=0.0035;Phy2=-0.003; Phy3 = 0.045;  
Pvy1=0.045;Pvy2=-0.03;Pvy3=-0.174;Pvy4=-0.45;
Rby1 = 6.38;

%----Lateral Force(pure side slip) input:横向侧偏 是侧偏角取tan后的值，即a*---%
r=0; %0.1*pi/180;                 %camber angle 侧倾角
sr=sin(r);                        %r 是外倾角，sr 表示r*
Cy=Pcy1;                          %取lamCy=1 
ta = tan(alpha *pi/180);          % 
uy=(Pdy1+Pdy2*dfz)*(1-Pdy3*sr^2);        %取lam(uy*)=1 
Dy=uy*Fz;   
Kya0=Pky1*Fz0*sin(2*atan(Fz/(Pky2*Fz0)));%取lam（Kya）=1 
Kya=Kya0*(1-Pky3*sr^2);
By=Kya/(Cy*Dy);                      %取伊布西隆y=0；  
Shy=Phy3*sr;
Svy=Fz*(Pvy3 + Pvy4*dfz)*sr;
ay=ta+Shy;  
Ey=(Pey1+Pey2*dfz)*(1-(Pey3+Pey4*sr)*sign(ay));
Fy0=Dy*sin(Cy*atan(By*ay-Ey.*(By*ay-atan(By*ay))))+Svy; 
% subplot(2,3,2) 
% plot(a,Fy0);grid  
% set(gca,'xlim',[-10 10])                      
% set(gca,'xtick',[-10:1:10]);                    
% set(gca,'ylim',[-4000 4000])                  
% set(gca,'ytick',[-4000:1000:4000]);             
% xlabel('侧偏角/（度）'); 
% ylabel('侧偏力/（N）'); 
% title('侧偏力(纯侧滑)');   


end % end of function

% for i= 1:1:50001
%     MF_Fy(i) = func_MFTyreModel(FzL1(i), alphaL1(i));
% end