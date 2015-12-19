function [reference_traj] = build_reference_traj(x_initial, x_final, N, restLengths, dt, testInitial)
% Build reference trajectory for spine based on desired final point x_f
% State Vector z refers to concatenation of state variabls x_k and input
% variables u_k: [x_1; x_2; ... ; x_N; u_1; u_2; ... ; u_N]

% Build initial trajectory to optimize around
z = x_initial;
z_f = x_final;
for k = 2:N
    x_k = ((k-1)/(N-1))*(x_final - x_initial) + x_initial;
    z = [z; x_k];
end

% Appending "reference" inputs to state vector z
for k = 1:(N-1)
    z = [z; (-1)^k*[0.01; -0.01; 0.01; -0.01; 0.01; -0.01; 0.01; -0.01]];
end

z = testInitial;

% Cost Function Definition
A = zeros(12);
A(1:6, 1:6) = 25*eye(6);
B = kron(eye(N), A);
Q = zeros(12*N + 8*(N-1));
Q(1:(12*N), 1:(12*N)) = B;
k = 12*(N-1)+1;
Q(k:(k+2), k:(k+2)) = 150*eye(3);

q = -Q*z;
q = q';

f = @(x) 0.5*z'*Q*z;

% Linear constraints
A_eq = zeros(12*N + 8*(N-1));
A_eq(1:12, 1:12) = eye(12);
b_eq = zeros(12*N + 8*(N-1), 1);
b_eq(1:12) = x_initial;

A_ineq = zeros(2*8*(N-1), 12*N + 8*(N-1));
A_ineq(1:(8*(N-1)), (12*N+1):end) = eye(8*(N-1));
A_ineq((8*(N-1)+1):end, (12*N+1):end) = -eye(8*(N-1));
b_ineq = 0.04*ones(2*8*(N-1), 1);

A_ineq1 = zeros(1, 12*N + 8*(N-1));
A_ineq1(1, (N-1)*12 + 6) = -1;
b_ineq1 = 0.02;

A_ineq = [A_ineq; A_ineq1];
b_ineq = [b_ineq; b_ineq1];
% Nonlinear constraints
g = @(x) 0;
h = @(x) trajectory_dynamics(x, N, restLengths, dt);

% Configuration parameters
user_cfg = struct();
user_cfg.min_approx_improve = 1e-1;
user_cfg.min_trust_box_size = 1e-3;
user_cfg.full_hessian = false;
user_cfg.h_use_numerical = false; % Provide numerical jach
user_cfg.initial_trust_box_size = 0.1;
user_cfg.max_merit_coeff_increases = 5;
user_cfg.initial_penalty_coeff = 1000;

reference_traj = penalty_sqp(z, Q, q, f, A_ineq, b_ineq, A_eq, b_eq, g, h, user_cfg);
end