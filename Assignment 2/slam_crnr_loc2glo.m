function [detected_corners_global] = slam_crnr_loc2glo(detected_corners, part_x, part_y, part_theta)
% Converts corners from local to global frame, preserves local x and y (legacy) and covariance

	if (~isempty(detected_corners))
		detected_corners_global = struct('x', 0, 'y', 0, 'heading', 0, 'angle', 0, 'local_x', 0, 'local_y', 0, 'covariance', zeros(4,4));
		for i = 1:size(detected_corners, 2)
			% Local to global transformation 
			detected_corners_global(i).x = detected_corners(i).x*cos(part_theta) - detected_corners(i).y*sin(part_theta) + part_x;
			detected_corners_global(i).y = detected_corners(i).x*sin(part_theta) + detected_corners(i).y*cos(part_theta) + part_y;
            
            detected_corners_global(i).heading = slam_in_pi(detected_corners(i).heading + part_theta);
            detected_corners_global(i).angle = detected_corners(i).angle;
            
			detected_corners_global(i).local_x = detected_corners(i).x;
			detected_corners_global(i).local_y = detected_corners(i).y;
            
			detected_corners_global(i).covariance = detected_corners(i).covariance;
            
            detected_corners_global(i).covariance(1, 1) = abs(detected_corners_global(i).covariance(1, 1) * cos(part_theta) - detected_corners_global(i).covariance(2, 2) * sin(part_theta)); 
            detected_corners_global(i).covariance(2, 2) = abs(detected_corners_global(i).covariance(1, 1) * sin(part_theta) + detected_corners_global(i).covariance(2, 2) * cos(part_theta)); 
		end
	else
		detected_corners_global = [];
	end
end