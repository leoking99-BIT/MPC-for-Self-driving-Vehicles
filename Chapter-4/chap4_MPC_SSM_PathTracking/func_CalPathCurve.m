function K=func_CalPathCurve(XA,YA,XB,YB,XC,YC)
%%
%求斜率
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

%% 最初的处理方法，逻辑不对，易出错。
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


function K=GetPathHeading(Xb,Yb,Xn,Yn)
    AngleY=Yn-Yb;
    AngleX=Xn-Xb;

%% 
    %***求Heading Angle 在0~2*pi之间 *******%
    if Xb==Xn
        if Yn>Yb
            K=pi/2;
        else
            K=3*pi/2;
        end
    else
        if Yb==Yn
            if Xn>Xb
                K=0;
            else
                K=pi;
            end
        else
            K=atan(AngleY/AngleX);
        end    
    end

%****修正K,使之在0~360°之间*****%
   if (AngleY>0&&AngleX>0)%第一象限
        K=K;
    elseif (AngleY>0&&AngleX<0)||(AngleY<0&&AngleX<0)%第二、三象限
        K=K+pi;
    else if (AngleY<0&&AngleX>0)%第四象限
            K=K+2*pi;  
        else
            K=K;
        end
   end
    
end % end of function

