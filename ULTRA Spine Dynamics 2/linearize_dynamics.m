function [A, B, c] = linearize_dynamics(xbar, ubar, restLengths, dt)
% Compute linearized dynamics based on linearization points xbar and ubar

eps = .001;
A = zeros(12);

for k = 1:12
    xlinU = xbar;
    xlinL = xbar;
    xlinU(k) = xbar(k) + eps;
    xlinL(k) = xbar(k) - eps;
    A(:, k) = (simulate_dynamics(xlinU, restLengths, ubar, dt, 0) - simulate_dynamics(xlinL, restLengths, ubar, dt, 0))/(2*eps);
end

B = zeros(12, 8);
for k = 1:8
    ulinU = ubar;
    ulinL = ubar;
    ulinU(k) = ubar(k) + eps;
    ulinL(k) = ubar(k) - eps;
    B(:, k) = (simulate_dynamics(xbar, restLengths, ulinU, dt, 0) - simulate_dynamics(xbar, restLengths, ulinL, dt, 0))/(2*eps);
end

c = simulate_dynamics(xbar, restLengths, ubar, dt, 0)' - A*xbar' - B*ubar;
end