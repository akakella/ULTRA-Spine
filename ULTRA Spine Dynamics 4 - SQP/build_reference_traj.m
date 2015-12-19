 function [reference_traj] = build_reference_traj(initial_traj, restLengths, dt, links, rad, obstacles)
% Build reference trajectory for spine based on desired final point x_f
% State Vector z refers to concatenation of state variabls x_k and input
% variables u_k: [x_1; x_2; ... ; x_N; u_1; u_2; ... ; u_N]

N = size(initial_traj, 2);
initial_traj = reshape(initial_traj(1:2:3, :), 2*N, 1);
%% Cost Function Definition
Q = 2*eye(2*N);

q = -Q*initial_traj;
q = q';

f = @(x) .5*initial_traj'*Q*initial_traj;

%% Linear constraints
A_eq = zeros(2*N);
A_eq(1:2, 1:2) = eye(2);
A_eq((end-1):end,(end-1):end) = eye(2);
b_eq = zeros(2*N, 1);
b_eq(1:2) = initial_traj(1:2);
b_eq((end-1):end) = initial_traj((end-1):end);

A_ineq = zeros(2*N);
b_ineq = zeros(2*N, 1);

%% Nonlinear constraints
% For enforcing obstacle constraints
make_robot_poly = @(x) orientedBoxToPolygon([x(1), x(2), rad, rad, 0]);
g = @(x) check_collisions(x, 0.02, [2, N], make_robot_poly, obstacles);

% For enforcing system dynamics
h = @(x) 0;
% h = @(x) trajectory_dynamics(x, N, restLengths, links, dt);

%% Configuration parameters
user_cfg = struct();
user_cfg.min_approx_improve = 1e-2;
user_cfg.min_trust_box_size = 1e-3;
user_cfg.full_hessian = false;
user_cfg.h_use_numerical = true; % Not providing numerical jach
user_cfg.initial_trust_box_size = .1;
user_cfg.max_merit_coeff_increases = 5;
user_cfg.initial_penalty_coeff = 1000;

reference_traj = penalty_sqp(initial_traj, Q, q, f, A_ineq, b_ineq, A_eq, b_eq, g, h, user_cfg);
end