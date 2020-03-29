function [y, Bout] = func_RLSEstimation_Ccr(f, d, FIR_Num, delta_n )
%----------------------------------------------------------%
% RLSFilt-Recursive Least-Squares FIR filter demonstration
% Usage : 
% 1) Initialization:
%     y = RLSFilt('initial', lambda_cr, Num_cr, delta)
%     d = Lambda: is the convergence rate parameter.
%     d = lambda_cr: is also called the "forgetting" exponential weight factor
%     Num_cr is the filter length
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
    persistent Y_cr F_cr B_cr lambda_cr  Num_cr  Rinv_cr  %delta
    % The following is initialization, and is executed once
    if (ischar(f) && strcmp(f,'initial')) % Initial
        lambda_cr = d;
%         delta = delta_n;
        Num_cr = FIR_Num;
        Rinv_cr = delta_n*eye(Num_cr);
        F_cr = zeros(Num_cr,1);
        Y_cr = zeros(1, Num_cr);
        B_cr = -92000; %zeros(Num_cr,1);
        y = 0;
        Bout = 0;
    else  % Filtering:
        for J = Num_cr:-1:2
            F_cr(J) = F_cr(J-1);
            Y_cr(J) = Y_cr(J-1);
        end;
        F_cr(1) = f;
        Y_cr(1) = d;
        % Perform the convolution
        y= F_cr'*B_cr;
        error=Y_cr-y;
        % Kalman gains
        K = Rinv_cr*F_cr/(lambda_cr + F_cr'*Rinv_cr*F_cr);
        % Update Rinv_cr
        Rinvn = (Rinv_cr - K*F_cr'*Rinv_cr)/lambda_cr;
        % Update the filter coefficients
        B_cr = B_cr + K * error;
        Rinv_cr = Rinvn;
        
        Bout = B_cr;
%         R_r_out = Rinv_cr;
    end
end