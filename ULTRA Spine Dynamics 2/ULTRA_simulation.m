clear variables
close all
clc
lim =0.475; %plot limits
dt = 0.01/2; %time step for dynamics
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
xbar = [0; 0; 0.1; 0; 0; 0; 0; 0; 0; 0; 0; 0];
ubar = [0; 0; 0; 0; 0; 0; 0; 0];

x_f = [0; 0; 0.3; 0; 0; 0; 0; 0; 0; 0; 0; 0]; % Desired State
e = 0;
%% Simulation Start
x=0; y=-0.0; z=0.0976; T=0.0; G=0.0; P=0;
dx=0; dy=0; dz=0; dT=0; dG=0; dP=0;

% 0.185
% 0.387
restLengths(1) = 0.1; % Vertical cable rest length
restLengths(2) = 0.1;
restLengths(3) = 0.1;
restLengths(4) = 0.1;
restLengths(5) = 0.187; % Saddle cable rest length
restLengths(6) = 0.187;
restLengths(7) = 0.187;
restLengths(8) = 0.187;

[A, B, c] = linearize_dynamics(xbar', ubar, restLengths, dt);

Figs = figure('units','normalized','outerposition',[0 0 1 1]);

M = struct('cdata', cell(1,round(length(time)/10)), 'colormap', cell(1,round(length(time)/10)));

leg = (l^2 - (h/2)^2)^.5;
Tetra2 = [(leg + rad/2), 0, -h/2; -(leg + rad/2), 0, -h/2; 0, (leg + rad/2), h/2; 0, -(leg + rad/2), h/2];
Tetra1 = Tetra2;
Tetra = [(leg + rad/2), 0, -h/2; -(leg + rad/2), 0, -h/2; 0, (leg + rad/2), h/2; 0, -(leg + rad/2), h/2];

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
systemStates = zeros(1,12);
flip = -1;

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
    
    systemStates = [x,y,z,T,G,P,dx,dy,dz,dT,dG,dP];
    %% Controller
    noise = 0;
    N = 3;

    clear inputs
    clear states
    inputs = sdpvar(8, N-1);
    states = sdpvar(12, N);

    constraints = [states(:, 1) == systemStates'];
    for k = 1:(N-2)
        constraints = [constraints, states(:, k+1) == A*states(:, k) + B*inputs(:, k) + c, ...
            -0.02*ones(8, 1) <= inputs(:, k) <= 0.02*ones(8, 1), ...
            norm(inputs(1:6, 1) - inputs(1:6, k), inf) <= 0.1];
    end
    constraints = [constraints, states(:, N) == A*states(:, N-1) + B*inputs(:, N-1) + c, -0.04*ones(8, 1) <= inputs(:, N-1) <= 0.04*ones(8, 1)];
    % constraints = [constraints, norm(states(1:6, 1) - states(1:6, N), inf) <= 0.1];

    objective = 0;
    for k = 1:(N-1)
        objective = objective + (15^k)*(states(5, k) - x_f(5))^2 + (10^k)*norm(inputs(:, k), 2);
    end
    objective = objective + (15^N)*(states(5, N) - x_f(5))^2;
    
     tic;
    diagnostics = optimize(constraints, objective, sdpsettings('solver', 'gurobi', 'verbose', 0));
    toc
    if (diagnostics.problem == 1)
        inputs = ubar;
    else
        states = value(states);
        inputs = value(inputs);
        inputs = inputs(:, 1)
    end
    
    approx = A*systemStates' + B*inputs + c;
    systemStates = simulate_dynamics(systemStates, restLengths, inputs, dt, noise);
    approx(6);
    systemStates(6);
    xbar = systemStates'
    ubar = inputs(:, 1)
    [A, B, c] = linearize_dynamics(xbar', ubar, restLengths, dt);
    toc
    input_arr(:, i) = inputs;
    angle(i) = systemStates(6);
    %% Controller End
    
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
end