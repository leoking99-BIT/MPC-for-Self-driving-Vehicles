function [y, Bout] = func_RLSEstimation_Ccf(f, d, FIR_Num, delta_n )
%----------------------------------------------------------%
% RLSFilt-Recursive Least-Squares FIR filter demonstration
% Usage : 
% 1) Initialization:
%     y = RLSFilt('initial', lambda_cf, Num_cf, delta)
%     d = Lambda: is the convergence rate parameter.
%     d = lambda_cf: is also called the "forgetting" exponential weight factor
%     Num_cf is the filter length
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
    persistent Y_cf F_cf B_cf lambda_cf  Num_cf  Rinv_cf  %delta
    % The following is initialization, and is executed once
    if (ischar(f) && strcmp(f, 'initial')) % Initial
        lambda_cf = d;
%         delta = delta_n;
        Num_cf = FIR_Num;
        Rinv_cf = delta_n*eye(Num_cf);
        F_cf = zeros(Num_cf,1);
        Y_cf = zeros(1, Num_cf);
        B_cf =  -92000; %zeros(1,1);
        y = 0;
        Bout = 0;
    else  % Filtering:
        for J = Num_cf:-1:2
            F_cf(J) = F_cf(J-1);
            Y_cf(J) = Y_cf(J-1);
        end;
        F_cf(1) = f;
        Y_cf(1) = d;
        % Perform the convolution
        y= F_cf'*B_cf;
        error=Y_cf-y;% 1*nd
        % Kalman gains
        K = Rinv_cf*F_cf/(lambda_cf + F_cf'*Rinv_cf*F_cf);
        % Update Rinv_cf
        Rinvn = (Rinv_cf - K*F_cf'*Rinv_cf)/lambda_cf;
        % Update the filter coefficients
        B_cf = B_cf + K* error;
        Rinv_cf = Rinvn;

        Bout = B_cf;     
%         R_f_out = Rinv_cf;
    end
end