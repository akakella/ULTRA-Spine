%% Initialize system parameters
clear variables
close all
clc

lim =0.675; % Plot limits
dt = 0.01; % Time step for dynamics
frame = 3; % Animation Frame Divisor
stringEnable = 1; % Enable string plotting
rad = 0.02; % Radius of Link
g = 9.81; % Gravity
h = 0.25; % Height of Link
l = 0.3; % Length of leg
leg = (l^2 - (h/2)^2)^.5; % Projection of leg length
links = 2; % Number of links in addition to base link

time = 0:dt:500; % Simulation time

%% Initialize Plot
Figs = figure('Units','Normalized');
M = struct('cdata', cell(1,round(length(time)/10)), 'colormap', cell(1,round(length(time)/10)));

cmaps = summer(512);
ax = axes();

grid on;
axis equal;
xlim([-lim-0.085,lim+0.085])
ylim([-lim,lim])
zlim([-0.15, lim])
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
colormap(cmaps(1:256,:))
shading interp
light
lighting phong
hold on
view(3)
title('Euler-Lagrange Model')
m = 256;

%% Simulation
restLengths(1) = 0.185; % Vertical cable rest length
restLengths(2) = 0.185;
restLengths(3) = 0.185;
restLengths(4) = 0.185;
restLengths(5) = 0.387; % Saddle cable rest length
restLengths(6) = 0.387;
restLengths(7) = 0.387;
restLengths(8) = 0.387;

systemStates = zeros(links, 12);

for k = 1:links
    x(k) = 0; y(k) = 0.0; z(k) = 0.18*k; T(k) = 0.0; G(k) = 0.0; P(k) = 0.0;
    dx(k) = 0; dy(k) = 0; dz(k) = 0; dT(k) = 0; dG(k) = 0; dP(k) = 0;
    
    Tetra{k} = [(leg + rad/2), 0, -h/2; -(leg + rad/2), 0, -h/2; 0, (leg + rad/2), h/2; 0, -(leg + rad/2), h/2];
    [transform{k}, handle{k}] = plotSpineLink(Tetra{k}, rad, ax);
end
Tetra{links+1} = Tetra{links};
[transform{links+1}, handle{links+1}] = plotSpineLink(Tetra{links+1}, rad, ax);

if (stringEnable)    
    anchor1=[0 0 rad];
    anchor2=[0 0 rad];
    String_pts = [];
    for k = 1:links
       String_pts = [String_pts; (Tetra{k}(1,:)+anchor1); (Tetra{k+1}(1,:)-anchor2); (Tetra{k}(3,:)+anchor1); (Tetra{k+1}(3,:)-anchor2); (Tetra{k}(3,:)+anchor1); (Tetra{k+1}(2,:)-anchor2); ...
                    (Tetra{k}(2,:)+anchor1); (Tetra{k+1}(2,:)-anchor2); (Tetra{k}(4,:)+anchor1); (Tetra{k+1}(4,:)-anchor2); (Tetra{k}(4,:)+anchor1); (Tetra{k+1}(1,:)-anchor2)];
    end
    string_handle=plot3(String_pts(:,1),String_pts(:,2),String_pts(:,3),'LineWidth',2,'Color','m');
    linkdata on
end

for (i = 2:length(time))
    String_pts = [];
    if mod(i, frame) == 0
        tic
        
        while toc < (dt*frame)
            % Wait until time has passed
        end
        
        for k = 1:links
            RR{k} =  getHG_Tform(x(k),y(k),z(k),T(k),G(k),P(k)); % Build graphical model of each link
            set(transform{k},'Matrix',RR{k});
        end
        
        if (stringEnable)
            for k = 2:(links+1)
                Tetra{k} = [(l^2 - (h/2)^2)^.5, 0, -h/2, 1; -(l^2 - (h/2)^2)^.5, 0, -h/2, 1; 0, (l^2 - (h/2)^2)^.5, h/2, 1; 0, -(l^2 - (h/2)^2)^.5, h/2, 1];
                Tetra{k} = RR{k-1}*Tetra{k}';
                Tetra{k} = Tetra{k}';
                Tetra{k} = Tetra{k}(:,1:3);
            end
            anchor2 = [0 0 rad];
            for k = 1:links
                String_pts = [String_pts; (Tetra{k}(1,:)+anchor1); (Tetra{k+1}(1,:)-anchor2); (Tetra{k}(3,:)+anchor1); (Tetra{k+1}(3,:)-anchor2); (Tetra{k}(3,:)+anchor1); (Tetra{k+1}(2,:)-anchor2); ...
                    (Tetra{k}(2,:)+anchor1); (Tetra{k+1}(2,:)-anchor2); (Tetra{k}(4,:)+anchor1); (Tetra{k+1}(4,:)-anchor2); (Tetra{k}(4,:)+anchor1); (Tetra{k+1}(1,:)-anchor2)];
            end
            refreshdata(string_handle);
        end
        drawnow;
    end
    
    for k = 1:links
        systemStates(k, 1) = x(k); systemStates(k, 2) = y(k); systemStates(k, 3) = z(k);
        systemStates(k, 4) = T(k); systemStates(k, 5) = G(k); systemStates(k, 6) = P(k);
        systemStates(k, 7) = dx(k); systemStates(k, 8) = dy(k); systemStates(k, 9) = dz(k);
        systemStates(k, 10) = dT(k); systemStates(k, 11) = dG(k); systemStates(k, 12) = dP(k);
    end
    
    noise = 1; % 1 for noise, 0 for no noise
    inputs = zeros(links, 8);
    systemStates = simulate_dynamics(systemStates, restLengths, inputs, dt, links, noise);
    
    for k = 1:links
        x(k) = systemStates(k, 1); y(k) = systemStates(k, 2); z(k) = systemStates(k, 3);
        T(k) = systemStates(k, 4); G(k) = systemStates(k, 5); P(k) = systemStates(k, 6);
        dx(k) = systemStates(k, 7); dy(k) = systemStates(k, 8); dz(k) = systemStates(k, 9);
        dT(k) = systemStates(k, 10); dG(k) = systemStates(k, 11); dP(k) = systemStates(k, 12);
    end
end
