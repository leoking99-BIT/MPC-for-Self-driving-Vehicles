 time = 1:1:2001;
 t = 0.05*time;
 num=length(t);
 Spet_speed = zeros(1,num);

 for index=1:num
    if index < 800
        Spet_speed(index) = 10;
    else
        if index < 1500
            Spet_speed(index) = 20;
        else
            Spet_speed(index) = 5;
        end
    end   
 end
 
 plot(time,Spet_speed)
 hold on
 grid on
 
 plot(u.signals.values(:,4),'r.');
