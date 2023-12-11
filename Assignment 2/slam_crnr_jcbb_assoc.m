function [assoc_corners, new_corners] = slam_crnr_jcbb_assoc(particle, detected_corners)
% Data association using Joint Compatibilty Branch and Bound
	
	assoc_corners = [];
	new_corners = [];
	detected_corners_count = size(detected_corners, 2);
    
	if (detected_corners_count ~= 0)
		if (particle.known_corners_count ~= 0)
			% Check for individual compatibility
			indiv_compat_matrix = zeros(detected_corners_count, particle.known_corners_count);
			for i = 1:detected_corners_count
				for j = 1:particle.known_corners_count

                    a = detected_corners(i).x - particle.corners(j).x;
					b = detected_corners(i).y - particle.corners(j).y;
					c = slam_in_pi(detected_corners(i).heading - particle.corners(j).heading);
					d = slam_in_pi(detected_corners(i).angle - particle.corners(j).angle);
                    
					Qt = detected_corners.covariance;
					covariance = particle.corners(j).covariance + Qt;
                    
					mahalanobis_dist = [a b c d] * (covariance^-1) * [a; b; c; d];
					if (mahalanobis_dist < particle.crnr_indiv_compat)
						indiv_compat_matrix(i,j) = 1;
					end
				end
			end
			
			% Start with hypothesis that all detected are new
			hypothesis = zeros(1, detected_corners_count);
			% Select candidates
			candidates_count = 1;
			candidates = struct('hypothesis', hypothesis, 'assoc_count', 0);
			while (~isempty(hypothesis))
				% Go to next hypothesis
				hypothesis = slam_hypothesis_next(hypothesis, indiv_compat_matrix);
				if (isempty(hypothesis))
					break;
				end
				% Count number of associations
				assoc_count = 0;
				for i = 1:detected_corners_count
					if (hypothesis(i) ~= 0)
						assoc_count = assoc_count + 1;
					end
				end
				% The more the merrier
				if (assoc_count > candidates(1).assoc_count)
					candidates_count = 1;
					candidates = struct('hypothesis', hypothesis, 'assoc_count', assoc_count);
				elseif (assoc_count == candidates(1).assoc_count)
					candidates_count = candidates_count + 1;
					candidates(candidates_count) = struct('hypothesis', hypothesis, 'assoc_count', assoc_count);
				end
			end
			
			% Decide on winner
			if (candidates_count > 1)
				% More than 1 candidate
				for i = 1:candidates_count
					joint_innovation = [];
					joint_covariance = [];
					for j = 1:detected_corners_count
						if (candidates(i).hypothesis(j) ~= 0)
							a = detected_corners(j).x - particle.corners(candidates(i).hypothesis(j)).x;
							b = detected_corners(j).y - particle.corners(candidates(i).hypothesis(j)).y;
							c = slam_in_pi(detected_corners(j).heading - particle.corners(candidates(i).hypothesis(j)).heading);  
							d = slam_in_pi(detected_corners(j).angle - particle.corners(candidates(i).hypothesis(j)).angle);
                            
							Qt = detected_corners(j).covariance;
							covariance = particle.corners(candidates(i).hypothesis(j)).covariance + Qt;
							joint_innovation = [joint_innovation a b c d];
							joint_covariance = blkdiag(joint_covariance, covariance);
						end
					end
					% Joint compatibility
					joint_mahalanobis_dist = joint_innovation*(joint_covariance^-1)*joint_innovation';
					if (i == 1)
						winner = candidates(1).hypothesis;
						winner_mahalanobis_dist = joint_mahalanobis_dist;
					else
						if (joint_mahalanobis_dist < winner_mahalanobis_dist)
							winner = candidates(i).hypothesis;
							winner_mahalanobis_dist = joint_mahalanobis_dist;
						end
					end
				end
			else
				winner = candidates(1).hypothesis;
			end
			% Break winning hypothesis into required outputs
			assoc_count = 0;
			new_count = 0;
			assoc_corners = struct('d_index', 0, 'k_index', 0);
			new_corners = struct('x', 0, 'y', 0, 'heading', 0, 'angle', 0, 'local_x', 0, 'local_y', 0, 'covariance', zeros(4,4));
			for i = 1:detected_corners_count
				if (winner(i) ~= 0)
					assoc_count = assoc_count + 1;
					assoc_corners(assoc_count).d_index = i;
					assoc_corners(assoc_count).k_index = winner(i);
				else
					new_count = new_count + 1;
					new_corners(new_count) = detected_corners(i);
				end
			end
			% In case no associations or no new corners
			if (assoc_count == 0)
				assoc_corners = [];
			end
			if (new_count == 0)
				new_corners = [];
			end
			% End of process
		else
			% No known corners, all detected are new
			new_corners = detected_corners;
		end
	else
		if (particle.known_corners_count ~= 0)
			% No corners detected
		end
	end
end

