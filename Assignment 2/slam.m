close all; clear all;
load Lidar_input;
load motion_estimate;

video_rec = 0; % Turn on/off video recording function 

%% SLAM Configuration
effective_range = 15000; % Use LiDAR range measurements shorter than this value (unit: mm)
particles_count = 50; % Total no. of particles
particle_init_weight = 10^100; % Initial weightage of each particle 
corners_indiv_compat = 300; % Individual compatibility threshold for corners
mem_blocks_size = 50; % For block memory allocation

translation_noise_add = 50; % Translation input noise (unit: mm)
theta_noise_add = 0.05; % Rotation input noise (unit: rad)

translation_noise_proportion = 0.05;
theta_noise_proportion = 0.05;

%% LiDAR sensor parameters
field_of_view = deg2rad(270); 
readings_count = 1081; 

%% Create auxiliary variable for later faster computation
bearing_mat = zeros(1, readings_count);
for i = 1:readings_count
	bearing_mat(i) = (pi-field_of_view)/2 + (i-1)*field_of_view/(readings_count-1);
	if (bearing_mat(i) > pi)
		bearing_mat(i) = bearing_mat(i)-2*pi;
	end
end
cos_bearing_mat = cos(bearing_mat);
sin_bearing_mat = sin(bearing_mat);

%% Create figure
figure_handle = figure('Name', 'SLAM', 'Position', get(0, 'Screensize'));
map_axes = subplot(1, 2, 2); % Where map is displayed
data_axes = subplot(1, 2, 1); % Where data is displayed

%% Prepare video recording
if (video_rec == 1)
    set(figure_handle, 'Renderer', 'zbuffer');
    video_out = VideoWriter('slam', 'MPEG-4');
    video_out.FrameRate = 10;
    open(video_out);
end

%% Initialise particles
for i = 1:particles_count
	known_corners(1:mem_blocks_size) = struct('x', 0, 'y', 0, 'heading', 0, 'angle', 0, 'covariance', zeros(4,4), 'view_count', 0);
    particles(i) = struct('x', 0, 'y', 0, 'theta', 0, 'weight', particle_init_weight, 'crnr_mem_blocks_count', 1, 'crnr_indiv_compat', corners_indiv_compat, 'known_corners_count', 0, 'corners', known_corners);
end

%% Limit lidar data range
output_temp = output;
for i = 1:size(output_temp,1)
	for j = 1:readings_count
		if (output_temp(i,j) > effective_range)
			output_temp(i,j) = 0;
		end
	end
end

%% Loop through data set
for i = 1:size(output_temp,1)
	%% Calculate Cartesian coordinates
	points = slam_points(output_temp(i,:), bearing_mat, cos_bearing_mat, sin_bearing_mat);
    
	%% Feature extraction, local frame
	[detected_lines_local detected_corners_local] = slam_lidar_feat_extrn(points);

    %% Motion prediction 
	for j = 1:particles_count		
        % Missing code start here ...
        
        rotation_noise_add = theta_noise_add * randn(1);
        rotation_noise_proportion = theta_noise_proportion * randn(1);
        delta_theta = (motion_estimate(j).theta * (1 + rotation_noise_proportion)) + rotation_noise_add;

        translation_noise_add_x = translation_noise_add * randn(1);
        translation_noise_proportion_x = translation_noise_proportion * randn(1);
        delta_x = (motion_estimate(j).x * (1 + translation_noise_proportion_x)) + translation_noise_add_x;

        translation_noise_add_y = translation_noise_add * randn(1);
        translation_noise_proportion_y = translation_noise_proportion * randn(1);
        delta_y = (motion_estimate(j).y * (1 + translation_noise_proportion_y)) + translation_noise_add_y;
        
        particles(j).theta = particles(j).theta + delta_theta;
        particles(j).x = particles(j).x + (cos(particles(j).theta)*delta_x) - (sin(particles(j).theta)*delta_y);
        particles(j).y = particles(j).y + (sin(particles(j).theta)*delta_x) + (cos(particles(j).theta)*delta_y);

        % Missing code end here ...
    end
    
    %% Update particle's map features and weight
	for j = 1:particles_count
		% Convert detected corners to global frame
		detected_corners_global = slam_crnr_loc2glo(detected_corners_local, particles(j).x, particles(j).y, particles(j).theta);
		% Corners association (in global frame)
		[assoc_corners new_corners] = slam_crnr_jcbb_assoc(particles(j), detected_corners_global);
		% Measurement update for associated corners 
		for k = 1:size(assoc_corners, 2)
			[particles(j).corners(assoc_corners(k).k_index) weight] = slam_crnr_kf(particles(j).corners(assoc_corners(k).k_index),detected_corners_global(assoc_corners(k).d_index));
			particles(j).weight = particles(j).weight * weight;
		end
		% Add new corners if not associated 
		particles(j) = slam_crnr_add(particles(j), new_corners);
    end 
    
    %% Draw result 
	slam_interface;
    
	%% Resample particles based on weights 
	particles = slam_resample(particles, particle_init_weight);
		
	%% Block memory allocation
	if (~mod(i, 5))
		for j = 1:particles_count
			if (particles(j).known_corners_count > 0.9*particles(j).crnr_mem_blocks_count*mem_blocks_size)
				particles(j).corners(particles(j).crnr_mem_blocks_count*mem_blocks_size + 1:(particles(j).crnr_mem_blocks_count+1)*mem_blocks_size) = struct('x', 0, 'y', 0, 'heading', 0, 'angle', 0, 'covariance', zeros(4, 4), 'view_count', 0);
				particles(j).crnr_mem_blocks_count = particles(j).crnr_mem_blocks_count + 1;
            end
		end
    end
    
end

%% Close video recording file 
if (video_rec == 1)
    close(video_out);
end