function [angle] = slam_in_pi(angle)
% Restrict heading angle value to be within [-pi, pi]
	
	while (angle > pi)
        angle = angle - 2*pi;
    end
    
    while (angle < -pi)
        angle = angle + 2*pi;
    end
    
end