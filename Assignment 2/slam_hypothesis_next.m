function [ next ] = slam_hypothesis_next( current, indiv_compat_matrix )
% Finds next valid hypothesis
	
	next = current;
	detected_count = size(next, 2);
	index = detected_count;
	
	known_count = size(indiv_compat_matrix, 2);
	repeats = 1;
	while (repeats == 1)
		% Do
		if (index ~= 0)
			next(index) = next(index) + 1;
			while (next(index) > known_count)
				next(index) = 0;
				index = index - 1;
				if (index == 0)
					break;
				end
				next(index) = next(index) + 1;
			end
		else
			next = [];
			break;
		end
		% Check again
		if (index == 0)
			next = [];
			break;
		end
		
		% While
		while (next(index) == 0 || (next(index) <= known_count && indiv_compat_matrix(index,next(index)) == 0))
			if (index ~= 0)
				next(index) = next(index) + 1;
				while (next(index) > known_count)
					next(index) = 0;
					index = index - 1;
					if (index == 0)
						break;
					end
					next(index) = next(index) + 1;
				end
			else
				next = [];
				break;
			end
			% Check again
			if (index == 0)
				next = [];
				break;
			end
		end
		% Check for repeats
		repeats = 0;
		if (~isempty(next) && detected_count >= 2)
			for i = 1:detected_count-1
				for j = i+1:detected_count
					if (next(i) ~= 0 && next(i) == next(j))
						repeats = 1;
						break;
					end
				end
				if (repeats == 1)
					break;
				end
			end
		end
	end
end