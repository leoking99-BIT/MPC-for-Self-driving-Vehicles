function [y, Bout, R_r_out] = func_RLSEstimation_Clr(f, d, FIR_Num, delta_n )
%----------------------------------------------------------%
% RLSFilt-Recursive Least-Squares FIR filter demonstration
% Usage : 
% 1) Initialization:
%     y = RLSFilt('initial', lambda_lr, Num_lr, delta)
%     d = Lambda: is the convergence rate parameter.
%     d = lambda_lr: is also called the "forgetting" exponential weight factor
%     Num_lr is the filter length
%     delta are the initial diagonal R^{-1}(n) matrix elements.
%     Example:
%     [y, e] = RLSFilt('initial', .95, 51, 0.01);
%     Note: RLSFilt returnsy=0 for initialization
% 2) Filtering:
%     [y, b] = RLSFilt(f, d);
%     where f is a single input value,
%     d is the desired value, and
%     y is the computed output value,
%     b is the coefficient vector.
%----------------------------------------------------------%
    persistent Y_lr F_lr B_lr lambda_lr  Num_lr  Rinv_lr  %delta
    % The following is initialization, and is executed once
    if (ischar(f) && strcmp(f,'initial')) % Initial
        lambda_lr = d;
%         delta = delta_n;
        Num_lr = FIR_Num;
        Rinv_lr = delta_n*eye(Num_lr);
        F_lr = zeros(Num_lr,1);
        Y_lr = zeros(1, Num_lr);
        B_lr = -92000; %zeros(Num_lr,1);
        y = 0;
        Bout = 0;
    else  % Filtering:
        for J = Num_lr:-1:2
            F_lr(J) = F_lr(J-1);
            Y_lr(J) = Y_lr(J-1);
        end;
        F_lr(1) = f;
        Y_lr(1) = d;
        % Perform the convolution
        y= F_lr'*B_lr;
        error=Y_lr-y;
        % Kalman gains
        K = Rinv_lr*F_lr/(lambda_lr + F_lr'*Rinv_lr*F_lr);
        % Update Rinv_lr
        Rinvn = (Rinv_lr - K*F_lr'*Rinv_lr)/lambda_lr;
        % Update the filter coefficients
        B_lr = B_lr + K * error;
        Rinv_lr = Rinvn;
        
        Bout = B_lr;
        R_r_out = Rinv_lr;
    end
end