 index = 1:1:2001;
 t = 0.05*index;
 
 A = 10;
 offst = 20; %正弦速度曲线的偏移
 T =50;
 f = 1/T;
 Ts = 0.05;
 y = A*sin(2*pi*f*t)+offst;
 plot(index,y)
 hold on
 
plot(u.signals.values(:,4),'r.');