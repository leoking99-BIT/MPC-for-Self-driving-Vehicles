function [y, Bout] = func_RLSFilter_Calpha_r_useFyr(f, d, FIR_Num, delta_n )
%----------------------------------------------------------%
% RLSFilt-Recursive Least-Squares FIR filter demonstration
% Usage : 
% 1) Initialization:
%     y = RLSFilt('initial', lambda_r, Num_r, delta)
%     d = lambda_r: is the convergence rate parameter,
%           also called the "forgetting" exponential weight factor
%     FIR_Num is the filter length
%     delta_n are the initial diagonal R^{-1}(n) matrix elements.
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
    persistent F_r B_r lambda_r  Num_r  Rinv_r  %delta
    % The following is initialization, and is executed once
    if (ischar(f) && strcmp(f,'initial')) % Initial
        lambda_r = d;
%         delta = delta_n;
        Num_r = FIR_Num;
        Rinv_r = delta_n*eye(Num_r);
        F_r = zeros(Num_r,1);
        B_r = zeros(Num_r,1);
        y = 0;
        Bout = 0;
    else  % Filtering:
        for J = Num_r:-1:2
            F_r(J) = F_r(J-1);
        end;
        F_r(1) = f;
        % Perform the convolution
        y= F_r'*B_r;
        error=d-y;
        % Kalman gains
        K = Rinv_r*F_r/(lambda_r + F_r'*Rinv_r*F_r);
        % Update Rinv_r
        Rinvn = (Rinv_r - K*F_r'*Rinv_r)/lambda_r;
        % Update the filter coefficients
        B_r = B_r + K*error;
        Bout = B_r;
        Rinv_r = Rinvn;
    end
end