function [controller] = genSTLController(x_initial, u_initial, restLengths, dt)

%% Create STL controller for Tensegrity Spine
[A, B, c] = linearize_dynamics(x_initial', u_initial, restLengths, dt);

C = eye(size(A, 2));
D = 0;

SG = STLC_lti(A, B, zeros(12,1), C, zeros(12,8), zeros(12, 1)); 
SG.x0= x_initial;

%% Controller Initialization
% Time
SG.time = 0:dt:1; % time for the dynamics
SG.ts=dt; % sampling time for controller
SG.L=3;  % horizon (# of steps)
SG.nb_stages=1; % repeats time

% Bounds
SG.u_ub(:) = 0.04;
SG.u_lb(:) = -0.04;
SG.nu = 8;
SG.nx = 12;

%% STL formula

SG.stl_list = {'ev_[0, 0.01] alw_[0, Inf] (abs(x3(t) - 0.1) > 0.01)'};

SG.min_rob = 0.01;    
SG.bigM = 1;

%% running stuff
fprintf('Computing controller...');
tic;
controller = get_controller(SG,'robust'); % interval
run_deterministic(SG, controller);
end