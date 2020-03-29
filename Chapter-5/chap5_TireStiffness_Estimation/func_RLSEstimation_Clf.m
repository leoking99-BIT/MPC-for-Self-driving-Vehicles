function [y, Bout] = func_RLSEstimation_Clf(f, d, FIR_Num, delta_n )
%----------------------------------------------------------%
% RLSFilt-Recursive Least-Squares FIR filter demonstration
% Usage : 
% 1) Initialization:
%     y = RLSFilt('initial', lambda_lf, Num_lf, delta)
%     d = Lambda: is the convergence rate parameter.
%     d = lambda_lf: is also called the "forgetting" exponential weight factor
%     Num_lf is the filter length
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
    persistent Y_lf F_lf B_lf lambda_lf  Num_lf  Rinv_lf  %delta
    % The following is initialization, and is executed once
    if (ischar(f) && strcmp(f, 'initial')) % Initial
        lambda_lf = d;
        Num_lf = FIR_Num;
        Rinv_lf = delta_n*eye(Num_lf);
        F_lf = zeros(Num_lf,1);
        Y_lf = zeros(1, Num_lf);
        B_lf =  -92000; %zeros(1,1);
        y = 0;
        Bout = 0;
    else  % Filtering:
        for J = Num_lf:-1:2
            F_lf(J) = F_lf(J-1);
            Y_lf(J) = Y_lf(J-1);
        end;
        F_lf(1) = f;
        Y_lf(1) = d;
        % Perform the convolution
        y= F_lf'*B_lf;
        error=Y_lf-y;% 1*nd
        % Kalman gains
        K = Rinv_lf*F_lf/(lambda_lf + F_lf'*Rinv_lf*F_lf);
        % Update Rinv_lf
        Rinvn = (Rinv_lf - K*F_lf'*Rinv_lf)/lambda_lf;
        % Update the filter coefficients
        B_lf = B_lf + K* error;
        Rinv_lf = Rinvn;

        Bout = B_lf;     
    end
end