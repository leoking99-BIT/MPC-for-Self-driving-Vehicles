function [PPx,PPy,de]=func_GetProjectPoint(Xb,Yb,Xn,Yn,Xc,Yc)
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