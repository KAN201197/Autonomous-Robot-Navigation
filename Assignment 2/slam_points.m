function [points] = slam_points(lidar_output, bearing_mat, cos_bearing_mat, sin_bearing_mat)
% Convert LiDAR range and angle measurements to Cartesian coordinates 

	readings_count = size(lidar_output, 2);
	points(1:readings_count) = struct('x',0,'y',0,'theta',0,'r',0);
    
	points_x = lidar_output.*cos_bearing_mat;
	points_y = lidar_output.*sin_bearing_mat;
    
	for i=1:readings_count
		points(i).x = points_x(i);
		points(i).y = points_y(i);
		points(i).theta = bearing_mat(i);
		points(i).r = lidar_output(i);
    end
    
end