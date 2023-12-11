function [particle] = slam_crnr_add(particle, new_corners)
% Initialises new fields in corners
	if (~isempty(new_corners))
		for i = 1:size(new_corners, 2)
			particle.known_corners_count = particle.known_corners_count + 1;
			particle.corners(particle.known_corners_count).x = new_corners(i).x;
			particle.corners(particle.known_corners_count).y = new_corners(i).y;
			particle.corners(particle.known_corners_count).heading = new_corners(i).heading;
			particle.corners(particle.known_corners_count).angle = new_corners(i).angle;
			particle.corners(particle.known_corners_count).covariance = new_corners(i).covariance;
			particle.corners(particle.known_corners_count).view_count = 1;
			% More new corners, smaller weight 
			particle.weight = particle.weight * 10^-50;
		end
	end
end