function [h, jach] = trajectory_dynamics(x, N, restLengths, dt)
nX = 12;
nU = 8;
x_traj = reshape(x(1:N*nX),nX, N);
u_traj = reshape(x(N*nX+1:end), nU, N-1);

h = zeros(nX*N + nU*(N-1), 1);
for i = 2:N
    k = nX*(i-1)+1;
    h(k:(k+nX-1)) = x_traj(:, i) - simulate_dynamics(x_traj(:, i-1)', restLengths, u_traj(:, i-1), dt, 0)';
end

jach = zeros(nX*N + nU*(N-1), nX*N + nU*(N-1));
for i = 2:N
    [A, B, ~] = linearize_dynamics(x_traj(:, i-1)', u_traj(:, i-1), restLengths, dt);
    
    f = nX*(i-1)+1; % Indexing for gradient of x(k+1) term
    g = nU*(i-1)+1;
    
    k = nX*(i-2)+1; % Indexing for gradient of f(x(k), u(k), dt) term
    m = nU*(i-2)+1;
    
    jach(k:(k+nX-1), k:(k+nX-1)) = -A;
    jach(k:(k+nX-1), (m+(nX*N)):(m+(nX*N)+nU-1)) = -B;
    
    jach(k:(k+nX-1), f:(f+nX-1)) = eye(nX);
end

end