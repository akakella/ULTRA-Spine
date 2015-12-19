function [refTraj] = referenceTrajectory(x0, xf)
    
    for k = 1:100
        x = x0;
        x(1) = x(1) + .001;
        x(2) = x(2) + .001;
        refTraj{k} = x;
    end




end