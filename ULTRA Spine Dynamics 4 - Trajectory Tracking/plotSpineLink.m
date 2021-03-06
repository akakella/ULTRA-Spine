function [transform,h] = plotSpineLink(spineCoord,rad,ax)
r = rad*ones(40,1);

% Plotting Center sphere of link
spineMean = mean(spineCoord, 1);
[x, y, z] = sphere(20);
x_c = rad*2*x + spineMean(1);
y_c = rad*2*y + spineMean(2);
z_c = rad*2*z + spineMean(3);

[x0, y0, z0] = cylinder2P(r,20,spineCoord(1,:),spineMean)
[x1, y1, z1] = cylinder2P(r,20,spineCoord(2,:),spineMean)
[x2, y2, z2] = cylinder2P(r,20,spineCoord(3,:),spineMean)
[x3, y3, z3] = cylinder2P(r,20,spineCoord(4,:),spineMean)


[x, y, z] = sphere;
x4 = rad*x + spineCoord(1,1); y4 = rad*y + spineCoord(1,2); z4 = rad*z + spineCoord(1,3);
x5 = rad*x + spineCoord(2,1); y5 = rad*y + spineCoord(2,2); z5 = rad*z + spineCoord(2,3);
x6 = rad*x + spineCoord(3,1); y6 = rad*y + spineCoord(3,2); z6 = rad*z + spineCoord(3,3);
x7 = rad*x + spineCoord(4,1); y7 = rad*y + spineCoord(4,2); z7 = rad*z + spineCoord(4,3);

hold on
h(1) = surf(ax,x_c,y_c,z_c,'LineStyle', 'none');
h(2) = surf(ax,x0,y0,z0,'LineStyle', 'none');
h(3) = surf(ax,x1,y1,z1,'LineStyle', 'none');
h(4) = surf(ax,x2,y2,z2,'LineStyle', 'none');
h(5) = surf(ax,x3,y3,z3,'LineStyle', 'none');
h(6) = surf(ax,x4,y4,z4,'LineStyle', 'none');
h(7) = surf(ax,x5,y5,z5,'LineStyle', 'none');
h(8) = surf(ax,x6,y6,z6,'LineStyle', 'none');
h(9) = surf(ax,x7,y7,z7,'LineStyle', 'none');
TT = hgtransform('Parent',ax);

transform = hgtransform('Parent',TT);
for i = 1:9
set(h(i),'Parent',transform)
end
end

%set(transform,'Matrix',Rflip*T2*Rx)

function [X, Y, Z] = cylinder2P(R, N,r1,r2)

    % The parametric surface will consist of a series of N-sided
    % polygons with successive radii given by the array R.
    % Z increases in equal sized steps from 0 to 1.

    % Set up an array of angles for the polygon.
    theta = linspace(0,2*pi,N);

    m = length(R);                 % Number of radius values
                                   % supplied.

    if m == 1                      % Only one radius value supplied.
        R = [R; R];                % Add a duplicate radius to make
        m = 2;                     % a cylinder.
    end


    X = zeros(m, N);             % Preallocate memory.
    Y = zeros(m, N);
    Z = zeros(m, N);
    
    v=(r2-r1)/sqrt((r2-r1)*(r2-r1)');    %Normalized vector;
    %cylinder axis described by: r(t)=r1+v*t for 0<t<1
    R2=rand(1,3);              %linear independent vector (of v)
    x2=v-R2/(R2*v');    %orthogonal vector to v
    x2=x2/sqrt(x2*x2');     %orthonormal vector to v
    x3=cross(v,x2);     %vector orthonormal to v and x2
    x3=x3/sqrt(x3*x3');
    
    r1x=r1(1);r1y=r1(2);r1z=r1(3);
    r2x=r2(1);r2y=r2(2);r2z=r2(3);
    vx=v(1);vy=v(2);vz=v(3);
    x2x=x2(1);x2y=x2(2);x2z=x2(3);
    x3x=x3(1);x3y=x3(2);x3z=x3(3);
    
    time=linspace(0,1,m);
    for j = 1 : m
      t=time(j);
      X(j, :) = r1x+(r2x-r1x)*t+R(j)*cos(theta)*x2x+R(j)*sin(theta)*x3x; 
      Y(j, :) = r1y+(r2y-r1y)*t+R(j)*cos(theta)*x2y+R(j)*sin(theta)*x3y; 
      Z(j, :) = r1z+(r2z-r1z)*t+R(j)*cos(theta)*x2z+R(j)*sin(theta)*x3z;
    end
end