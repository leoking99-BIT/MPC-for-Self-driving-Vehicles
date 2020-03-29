%***************************************************************%
% 这里生产双移线试验的参考道路
% 设定纵向速度为1m/s, 采样时间为0.1s, 因此每个路点在x方向上的距离为0.1m
% 产生的参考路径保存的.mat文件，[X_ref, Yref, Heading_ref]
%---------------------------------------------------------------%
% Published by: Kai Liu
% Email:leoking1025@gmail.com
% My homepage: https://sites.google.com/site/kailiumiracle/ 
%***************************************************************%
    %---------------双移线轨迹形状参数-------%
    shape = 2.4;%参数名称，用于参考轨迹生成
    dx1 = 25; 
    dx2 = 21.95;%没有任何实际意义，只是参数名称
    dy1 = 4.05; 
    dy2 = 5.7;%没有任何实际意义，只是参数名称
    Xs1 = 27.19; 
    Xs2 = 56.46;%参数名称

    DataNum = 3000; %
    DLS_path_cell = cell(DataNum,1);
    Ts    = 0.1;   
    X_DOT = 1.0;  %惯性坐标系下纵向速度
    
%%    
    X_0   = -50;  % 参考路径在X轴的起始坐标
    Line_segment_Num = 500;
    for p = 1 : 1 : Line_segment_Num
        X_ref       = X_0 + X_DOT * p * Ts; %首先计算出未来X的位置
        Y_ref       = 0;
        Heading_ref = 0;
        DLS_path_cell{p,1} = [X_ref, Y_ref, Heading_ref];
    end

    X_0 = 0;
    for p = 1 : 1 : DataNum-Line_segment_Num
        X_ref       = X_0 + X_DOT * p * Ts; %首先计算出未来X的位置
        z1          = shape/dx1*(X_ref - Xs1) - shape/2;
        z2          = shape/dx2*(X_ref - Xs2) - shape/2;
        Y_ref       = dy1/2*(1+tanh(z1)) - dy2/2*(1+tanh(z2));
        Heading_ref = atan(dy1*(1/cosh(z1))^2*(1.2/dx1) - dy2*(1/cosh(z2))^2*(1.2/dx2));
        DLS_path_cell{Line_segment_Num + p,1} = [X_ref, Y_ref, Heading_ref];
    end
    
    
    DLS_path=cell2mat(DLS_path_cell);

    save Waypoints_Double_Line_Shift.mat DLS_path;

%%
figure(1)
plot(DLS_path(:,1), DLS_path(:,2), 'k');

figure(2)
plot(1:DataNum, DLS_path(:,3), 'b');







