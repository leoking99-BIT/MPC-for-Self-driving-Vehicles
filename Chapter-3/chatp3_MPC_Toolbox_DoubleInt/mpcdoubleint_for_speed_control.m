%% Control of a Single-Input-Single-Output Plant
% This example shows how to control a double integrator plant under input
% saturation in Simulink(R).

% Copyright 1990-2012 The MathWorks, Inc.
% $Revision: 1.1.4.12 $  $Date: 2012/12/20 17:33:53 $   

%% MPC Controller Setup 
% Create MPC controller in the workspace.
Ts = .1;                                    % Sampling time
p = 20;                                     % Prediction horizon
m = 3;                                      % Control horizon
% mpc_controller = mpc(tf(1,[1 0 0]),Ts,p,m); % MPC object
% x_dot = Ax+Bu
% y = Cx +Du
A=[0 1; 0 0];
B=[0; 1];
C = [1 0];
D = [0];
mpc_controller = mpc(ss(A,B,C,D),Ts,p,m); % MPC object
mpc_controller.MV=struct('Min',-1,'Max',1); % Input saturation constraints

%% MPC Simulation Using Simulink(R)
% The example uses Simulink(R).
if ~mpcchecktoolboxinstalled('simulink')
    disp('Simulink(R) is required to run this example.')
    return
end

%%
% Setup simulation parameters.
x01=0;                                      % Initial state: First integrator
x02=0;                                      % Initial state: Second integrator
Tstop=5;                                    % Simulation time
r=2;                                        % Set point

%% 
% Run simulation.
open_system('mpc_doubleint_for_speedctrl');               % Open Simulink(R) Model
sim('mpc_doubleint_for_speedctrl',Tstop);                 % Start Simulation
    
