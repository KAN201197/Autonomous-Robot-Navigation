function [particles] = slam_resample(particles, init_weight)
	
	particles_count = size(particles, 2);
	weight_total = 0;
	
    for i = 1:particles_count
		weight_total = weight_total + particles(i).weight;
    end
    
     % Create a new array to store the resampled particles
    particles_new(particles_count) = struct('x', [], 'y', [], 'theta', [], 'weight', [], 'crnr_mem_blocks_count', [], 'crnr_indiv_compat', [], 'known_corners_count', [], 'corners', []);

	for i = 1:particles_count
        % Missing codes start here
        
        % Resamples particles based on their weights

        % Generate a random number between 0 and the total weight
        W_random = rand() * weight_total;

        W = weight_total;

        for j = 1:particles_count
            W = W - particles(j).weight;
            if W <= W_random

                % Particle i in the old set is chosen as one element in the new set
                particles_new(i) = particles(j);
                break;
            end
        end

        % Afterwards, each new partical should be given the same init_weight
        particles_new(i).weight = init_weight;

        % Missing codes end here
    end

    % Update particles with the new set
    particles = particles_new;

end