function [val,jac] = check_collisions(x, dsafe, traj_shape, make_robot_poly, obstacles)
assert(size(x,2)==1);
traj = reshape(x, traj_shape);
[K,T] = size(traj);

val = zeros(length(obstacles)*size(traj,2),1);
jac = zeros(size(val,1), size(x,1));

icontact = 1;


for t=1:T
    xt = traj(:,t);
    for iobs=1:length(obstacles)
        [d,pts] = signedDistancePolygons(...
                make_robot_poly(xt), ...
                obstacles{iobs});
        ptOnRobot = pts(1,:);
        ptOnObs = pts(2,:);

        gradd = zeros(1, K);
        gradd(1) = (signedDistancePolygons(make_robot_poly([xt(1)+.001; xt(2)]), obstacles{iobs}) - signedDistancePolygons(make_robot_poly([xt(1)-.001; xt(2)]), obstacles{iobs}))/.002;
        gradd(2) = (signedDistancePolygons(make_robot_poly([xt(1); xt(2)+.001]), obstacles{iobs}) - signedDistancePolygons(make_robot_poly([xt(1); xt(2)-.001]), obstacles{iobs}))/.002;
        val(icontact) = dsafe - d;
        jac(icontact,K*(t-1)+1:K*t) = gradd;
        icontact = icontact+1; 
   end 
end


end