% Bezier interpolation for given four control points.
% Each control point can be in N-Dimensional vector space.
% Input:
% P0,P1,P2,P3: four control points of bezier curve,
%              control points can have any number of coordinates
% t(optional arg):vector that holds paramter t values b/w 0 and 1 at which 
%                 bezier curve is evaluated (default 101 values between 0
%                 and 1.)

% Output:
% Q evaluated values of bezier curves. Number of columns of Q are equal to
% number of coordinates in control point. For example for 2-D, Q has two
% columns. Column 1 for x value and column 2 for y values. Similarly for
% 3-D, Q will have three columns

function Q=func_bezierInterp(P0,P1,P2,P3,varargin)

%%% Default Values %%%
t=linspace(0,1,101); % uniform parameterization 
defaultValues = {t};
%%% Assign Valus %%%
nonemptyIdx = ~cellfun('isempty',varargin);
defaultValues(nonemptyIdx) = varargin(nonemptyIdx);
[t] = deal(defaultValues{:});
% % --------------------------------
M=[-1  3 -3 1;
    3 -6  3 0;
   -3  3  0 0;
    1  0  0 0];
for k=1:length(t)
    Q(k,:)=[t(k)^3 t(k)^2 t(k) 1]*M*[P0;P1;P2;P3];
end
% % Ref: Mathematical Elements of Computer Graphics by
% %      David F. Rogers and J. Alan Adams (pg. 296)
% % --------------------------------
% % OR
% % Equation of Bezier Curve, utilizes Horner's rule for efficient computation.
% % Q(t)=(-P0 + 3*(P1-P2) + P3)*t^3 + 3*(P0-2*P1+P2)*t^2 + 3*(P1-P0)*t + Px0
% c3 = -P0 + 3*(P1-P2) + P3;
% c2 = 3*(P0 - (2*P1)+P2); 
% c1 = 3*(P1 - P0);
% c0 = P0;
% for k=1:length(t)
%     Q(k,:)=((c3*t(k)+c2)*t(k)+c1)*t(k) + c0;    
% end

% % % --------------------------------
% % % Author: Dr. Murtaza Khan
% % % Email : drkhanmurtaza@gmail.com
% % % --------------------------------
