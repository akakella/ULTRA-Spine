clear variables
close all
clc
lim = 0.275; %plot limits
dt = 0.005; %time step for dynamics
frames = 3; %frame divisor
stringEnable =1; %enable string plotting
rad = 0.01; % Radius of Link
g = 9.81; % Gravity
h = 0.25; % Height of Link
l = 0.3; % Length of spine "leg"

h = 0.15;
l = 0.15;
m = (0.282/5) * 1.2;

time = 0:dt:500; % Simulation time
xbar = [0; 0; 0.0976; 0; 0; 0; 0; 0; 0; 0; 0; 0];
ubar = [0; 0; 0; 0; 0; 0; 0; 0];

x_final = [0.0; 0.0; 0.07; 0; 0; 0.25; 0; 0; 0; 0; 0; 0]; % Desired State

%% Initial Conditions
x=0; y=-0.0; z=0.0976; T=0.0; G=0.0; P=0;
dx=0; dy=0; dz=0; dT=0; dG=0; dP=0;

restLengths(1) = 0.1; % Vertical cable rest length
restLengths(2) = 0.1;
restLengths(3) = 0.1;
restLengths(4) = 0.1;
restLengths(5) = 0.187; % Saddle cable rest length
restLengths(6) = 0.187;
restLengths(7) = 0.187;
restLengths(8) = 0.187;

% reference_traj = build_reference_traj(xbar, x_final, 50, restLengths, dt);
% reference_traj = reshape(reference_traj(1:600), 12, 50);

Figs = figure('units','normalized','outerposition',[0 0 1 1]);

M = struct('cdata', cell(1,round(length(time)/10)), 'colormap', cell(1,round(length(time)/10)));

[A, B, c] = linearize_dynamics(xbar', ubar, restLengths, dt);

leg = (l^2 - (h/2)^2)^.5;
Tetra2 = [(leg + rad/2), 0, -h/2; -(leg + rad/2), 0, -h/2; 0, (leg + rad/2), h/2; 0, -(leg + rad/2), h/2];
Tetra1 = Tetra2;
Tetra = [(leg + rad/2), 0, -h/2; -(leg + rad/2), 0, -h/2; 0, (leg + rad/2), h/2; 0, -(leg + rad/2), h/2];

%% Controller Initialization
N = 3;
inputs = sdpvar(repmat(8, 1, N-1), repmat(1, 1, N-1));
states = sdpvar(repmat(12, 1, N), repmat(1, 1, N));
A_t = sdpvar(repmat(12, 1, 12), repmat(1, 1, 12));
B_t = sdpvar(repmat(12, 1, 8), repmat(1, 1, 8));
c_t = sdpvar(12, 1);

input_lim = 0.05*ones(8, 1);
constraints = [];

for k = 1:(N-2)
    constraints = [constraints, states{k+1} == [A_t{:}]*states{k} + [B_t{:}]*inputs{k} + c_t, ...
        -input_lim <= inputs{k} <= input_lim, ...
        norm(inputs{k}(1:6) - inputs{1}(1:6)) <= 0.1];
end
constraints = [constraints, states{N} == [A_t{:}]*states{N-1} + [B_t{:}]*inputs{N-1} + c_t, -input_lim <= inputs{N-1} <= input_lim];
constraints = [constraints, norm(states{1}(1:6) - states{N}(1:6), inf) <= 0.1];

objective = 15*norm(states{1}(6) - x_final(6), 2);
for k = 2:(N-1)
    objective = objective + (15^k)*norm(states{k}(6) - x_final(6), 2) + (2^k)*norm(inputs{k}, 2) + (4^k)*norm(inputs{k} - inputs{k-1}, 2);
end
objective = objective + (15^N)*norm(states{N}(6) - x_final(6), 2);

parameters_in = {states{1}, [A_t{:}], [B_t{:}], c_t};
solutions_out = {[inputs{:}], [states{:}]};

controller = optimizer(constraints, objective, sdpsettings('solver', 'gurobi'), parameters_in, solutions_out);

%% Simulation
cmaps = summer(512);
ax = axes();
[T_base,h_base] = plotSpineLink(Tetra,rad,ax);
[T_top,h_top] = plotSpineLink(Tetra,rad,ax);
if(stringEnable)
    
    anchor1=[0 0 rad];
    anchor2=[0 0 rad];
    
    String_pts = [...
            (Tetra1(1,:)+anchor1);
            (Tetra2(1,:)-anchor2);
            (Tetra1(3,:)+anchor1);
            (Tetra2(3,:)-anchor2);
            (Tetra1(3,:)+anchor1);
            (Tetra2(2,:)-anchor2);
            (Tetra1(2,:)+anchor1);
            (Tetra2(2,:)-anchor2);
            (Tetra1(4,:)+anchor1);
            (Tetra2(4,:)-anchor2);
            (Tetra1(4,:)+anchor1);
            (Tetra2(1,:)-anchor2)];
    string_handle=plot3(String_pts(:,1),String_pts(:,2),String_pts(:,3),'LineWidth',2,'Color','m');
    linkdata on
end
grid on;
axis equal;
xlim([-lim-0.085,lim+0.085])
ylim([-lim,lim])
zlim([-0.15,lim])
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
colormap(cmaps(1:256,:))
shading interp
light
lighting phong
hold on
view(3)
title('ULTRA Spine Model')
m = 256;  % Number of colors in the current

systemStates = [x,y,z,T,G,P,dx,dy,dz,dT,dG,dP];
xhat = [0; 0; 0.1; 0; 0; 0; 0; 0; 0; 0; 0; 0];
Q = .001*eye(12);
% R = .001*eye(6);
R = .001;
Vfilt = Q;
% C = zeros(6, 12);
% C(6, 6) = 1;
C = zeros(1, 12);
C(6) = 1;

for i = 2:length(time)
    
    if mod(i,frames) == 0
        tic
        while toc < dt*frames
        end
        
        RR1 =  getHG_Tform(x,y,z,T,G,P);
        set(T_top,'Matrix',RR1);
        
        if stringEnable
            anchor2 = [0 0 rad 1];
            Tetra2 = [(l^2 - (h/2)^2)^.5, 0, -h/2, 1; -(l^2 - (h/2)^2)^.5, 0, -h/2, 1; 0, (l^2 - (h/2)^2)^.5, h/2, 1; 0, -(l^2 - (h/2)^2)^.5, h/2, 1];
            Tetra2 = RR1*Tetra2';
            Tetra2 = Tetra2';
            RR1(1:3,4) = [0; 0; 0];
            anchor2 = RR1*anchor2';
            anchor2 = anchor2(1:3)';
            Tetra2 = Tetra2(:,1:3);
            String_pts = [...
            (Tetra1(1,:)+anchor1);
            (Tetra2(1,:)-anchor2);
            (Tetra1(3,:)+anchor1);
            (Tetra2(3,:)-anchor2);
            (Tetra1(3,:)+anchor1);
            (Tetra2(2,:)-anchor2);
            (Tetra1(2,:)+anchor1);
            (Tetra2(2,:)-anchor2);
            (Tetra1(4,:)+anchor1);
            (Tetra2(4,:)-anchor2);
            (Tetra1(4,:)+anchor1);
            (Tetra2(1,:)-anchor2)];
            refreshdata(string_handle);
        end
        drawnow;
    end
    
    %% Controller
    tic;
    noise = 1;
    outputs = controller{{xhat, A, B, c}}; 
    control = outputs{1}(:, 1);
    systemStates = simulate_dynamics(systemStates, restLengths, control, dt, noise);
    %% State Update
    
    x = systemStates(1);
    y = systemStates(2);
    z = systemStates(3);
    T = systemStates(4); % about x
    G = systemStates(5); % about y
    P = systemStates(6); % about z
    dx = systemStates(7);
    dy = systemStates(8);
    dz = systemStates(9);
    dT = systemStates(10);
    dG = systemStates(11);
    dP = systemStates(12);
    
    systemStates = [x,y,z,T,G,P,dx,dy,dz,dT,dG,dP];
    
    xpred = simulate_dynamics(xhat', restLengths, control, dt, noise);
    Vk =  A*Vfilt*A' + Q; % V(t+1)|0:t
    
    K = Vk*C'*inv(C*Vk*C' + R);
    xhat = xpred' + K*(systemStates(6)' - C*xpred'); % x(t+1)|0:(t+1)
    Vfilt = (eye(12) - K*C)*Vk; % V(t+1)|0:(t+1)
    
    [A, B, c] = linearize_dynamics(xhat', control, restLengths, dt);
    
    truePosition(i, 1:6) = systemStates(1:6);
    estPosition(i, 1:6) = xhat(1:6);
    toc
end