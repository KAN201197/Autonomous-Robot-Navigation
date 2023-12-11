function [result_list] = slam_lidar_split_merge(points, threshold)
% Split and merge algorithm for line extraction
	
	last = size(points, 2); % Total no. of points
	
    x1 = points(1).x;
    y1 = points(1).y;
    x2 = points(last).x;
    y2 = points(last).y;
    
	if (abs(x2-x1) > abs(y2-y1))
		% y = mx + c
		m = (y2-y1) / (x2-x1); % Gradient
		c = y1 - m*x1; % y-intercept

        if (c >= 0)
            theta = atan(m) + pi/2;
        else
            theta = atan(m) - pi/2;
        end
        
        r = abs(c*cos(theta));
        theta = slam_in_pi(theta);

	else
		% x = my + c
		m = (x2-x1)/(y2-y1); % Gradient
		c = x1 - m*y1; % x-intercept

        if (c >= 0)
            theta = atan(m) + pi/2;
        else
            theta = atan(m) - pi/2;
        end
        
        r = abs(c*cos(theta));
        theta = slam_in_pi(pi/2 - theta);
        
    end
    
    % Missing codes start here ...
    
    % Find furthest point to this line 
	point_distances = zeros(1, last);
    winner_value = abs(point_distances(1));
    winner_index = 1;

    for i = 1:last

        % Calculate the distance of the point from the line
        distance = abs(((x2-x1)*(y1-points(i).y) - (x1-points(i).x)*(y2-y1))/sqrt((x2-x1)^2 + (y2-y1)^2));

        % store value into point_distance array
        point_distances(i) = distance;
        
        % comparison to find the highest distance of point to the line
        if distance >= winner_value
            winner_value = distance;
            winner_index = i;
        end
    end
      
    % Missing codes end here ...
	
	if (winner_value > threshold)
		% Split
		res_a = slam_lidar_split_merge(points(1: winner_index), threshold);
		res_b = slam_lidar_split_merge(points(winner_index: last), threshold);
		result_list = [res_a res_b];
	else
		%% Calculate / design covariance matrix for line
		r_sigma = sum(abs(point_distances))/last;
		line_covariance = [0.05 0.7; 0.7 r_sigma^2];
		
		l = sqrt((points(last).x-points(1).x)^2 + (points(last).y-points(1).y)^2);
		result_list = struct('p1', points(1), 'p2', points(last), 'theta', theta, 'r', r, 'length', l, 'covariance', line_covariance);
	end
end