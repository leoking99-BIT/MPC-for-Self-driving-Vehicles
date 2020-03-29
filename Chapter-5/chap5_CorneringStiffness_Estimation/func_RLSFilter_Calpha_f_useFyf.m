function [y, Bout] = func_RLSFilter_Calpha_f_useFyf(f, d, FIR_Num, delta_n )
%----------------------------------------------------------%
% RLSFilt-Recursive Least-Squares FIR filter demonstration
% Usage : 
% 1) Initialization:
%     y = RLSFilt('initial', lambda_f, Num_f, delta)
%     d = Lambda: is the convergence rate parameter.
%     d = lambda_f: is also called the "forgetting" exponential weight factor
%     Num_f is the filter length
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
    persistent Y_f F_f B_f lambda_f  Num_f  Rinv_f  %delta
    % The following is initialization, and is executed once
    if (ischar(f) && strcmp(f, 'initial')) % Initial
        lambda_f = d;
%         delta = delta_n;
        Num_f = FIR_Num;
        Rinv_f = delta_n*eye(Num_f);
        F_f = zeros(Num_f,1);
        Y_f = zeros(1, Num_f);
        B_f = -90000; % zeros(Num_f,1);
        y = 0;
        Bout = 0;
    else  % Filtering:
        for J = Num_f:-1:2
            F_f(J) = F_f(J-1);
            Y_f(J) = Y_f(J-1);
        end;
        F_f(1) = f;
        Y_f(1) = d;
        % Perform the convolution
        y= F_f'*B_f;
        error=Y_f-y; % 1*nd
        % Kalman gains
        K = Rinv_f*F_f/(lambda_f + F_f'*Rinv_f*F_f);
        % Update Rinv_f
        Rinvn = (Rinv_f - K*F_f'*Rinv_f)/lambda_f;
        % Update the filter coefficients
        B_f = B_f +error *K;
        Bout = B_f;
        Rinv_f = Rinvn;
    end
end