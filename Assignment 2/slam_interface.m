cla(data_axes); cla(map_axes);
set(data_axes, 'position', [0.04 0.1 0.5-3*0.02 0.8]);
set(map_axes, 'position', [0.5+0.02 0.1 0.5-3*0.02 0.8]);
hold on;

%% Data plots
subplot(data_axes);

% Draw angle coverage lines and range circles  
plot([0 effective_range*0.7087], [0 effective_range*-0.7056], 'k');
plot([0 -effective_range*0.7087], [0 effective_range*-0.7056], 'k');
for dist = 5000:5000:effective_range
	rectangle('Position',[-dist -dist dist*2 dist*2], 'Curvature', [1, 1]); 
end

% Draw data points
data_x = output_temp(i, :) .* cos_bearing_mat;
data_y = output_temp(i, :) .* sin_bearing_mat;
plot(data_x, data_y, 'r.');

% Draw local lines 
for j = 1:size(detected_lines_local, 2)
	plot([detected_lines_local(j).p1.x detected_lines_local(j).p2.x], [detected_lines_local(j).p1.y detected_lines_local(j).p2.y], 'b', 'LineWidth', 2);
end

% Draw local corners
for j = 1:size(detected_corners_local, 2)
	rectangle('Position', [detected_corners_local(j).x-250 detected_corners_local(j).y-250 500 500], 'Curvature', [1, 1]); % Circle
	text(detected_corners_local(j).x+400, detected_corners_local(j).y, num2str(rad2deg(detected_corners_local(j).angle), '%.1f'));
	plot([detected_corners_local(j).x detected_corners_local(j).x+1000*cos(detected_corners_local(j).heading)], [detected_corners_local(j).y detected_corners_local(j).y+1000*sin(detected_corners_local(j).heading)], 'g', 'LineWidth', 2);
end
axis equal;

%% Map plots
subplot(map_axes);

% Draw all particles and find the particle with the highest weight
weight_max = 0;
idx_max = 1;

for j = 1:particles_count
	plot(particles(j).x, particles(j).y, 'bx', 'MarkerSize', 10);
	plot([particles(j).x particles(j).x+1000*cos(particles(j).theta+0.5*pi)], [particles(j).y particles(j).y+1000*sin(particles(j).theta+0.5*pi)], 'b');
    
    if (particles(j).weight > weight_max)
        weight_max = particles(j).weight;
        idx_max = j;
    end
end

% Draw the partical with the highest weight 
plot(particles(idx_max).x, particles(idx_max).y, 'rx', 'MarkerSize',10);
plot([particles(idx_max).x particles(idx_max).x+1000*cos(particles(idx_max).theta+0.5*pi)], [particles(idx_max).y particles(idx_max).y+1000*sin(particles(idx_max).theta+0.5*pi)], 'r');

% Draw the highest weight particle's map features
for j = 1:particles(idx_max).known_corners_count
	% Corner depiction
	plot([particles(idx_max).corners(j).x particles(idx_max).corners(j).x-500*cos(particles(idx_max).corners(j).heading+0.5*particles(idx_max).corners(j).angle)],[particles(idx_max).corners(j).y particles(idx_max).corners(j).y-500*sin(particles(idx_max).corners(j).heading+0.5*particles(idx_max).corners(j).angle)],'b');
	plot([particles(idx_max).corners(j).x particles(idx_max).corners(j).x-500*cos(particles(idx_max).corners(j).heading-0.5*particles(idx_max).corners(j).angle)],[particles(idx_max).corners(j).y particles(idx_max).corners(j).y-500*sin(particles(idx_max).corners(j).heading-0.5*particles(idx_max).corners(j).angle)],'b');
	% Corner uncertainty depiction
    size_x = 10 * sqrt(particles(idx_max).corners(j).covariance(1, 1));
    size_y = 10 * sqrt(particles(idx_max).corners(j).covariance(2, 2));
	rectangle('Position', [particles(idx_max).corners(j).x-0.5*size_x particles(idx_max).corners(j).y-0.5*size_y size_x size_y], 'Curvature', [1, 1]);
end

%% 
axis equal;
hold off;
drawnow;

%% Write to video recording file 
if (video_rec == 1)
    writeVideo(video_out, getframe(figure_handle));
end