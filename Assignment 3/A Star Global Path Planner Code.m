% Load the map matrix (map data)
load('Map.mat');

% Define start and goal pose
% start = [1, 1];
% goal = [10, 10];

% Test generate a random start and goal pose
start = [randi(10), randi(10)];
goal = [randi(10), randi(10)];

% Check if the start pose is in an occupied area or not
if Map(start(1), start(2)) == 1
    disp('Start pose is in an occupied area. Cannot generate a path !');
else
    % Generate a path (waypoint) based on A* algorithm
    path = astarpathplanner(Map, start, goal);

    % Create a figure to visualize the map
    figure;

    % Display the binary map and static obstacles on the figure (white color
    % indicates free space, and black color indicates obstacles)
    imagesc(Map, 'CDataMapping', 'scaled');
    colormap([1, 1, 1; 0, 0, 0]);
    axis equal;
    title('Binary Map and Static Obstacles');

    hold on;

    % If a path is found, plot the path on the map
    if ~isempty(path)
        % Mark the start pose with a green star
        plot(start(2), start(1), 'g*', 'MarkerSize', 10);
        % Mark the goal pose with a blue star
        plot(goal(2), goal(1), 'b*', 'MarkerSize', 10);
        % Plot the path according to path data
        plot(path(:, 2), path(:, 1), 'r-o', 'LineWidth', 2, 'MarkerSize', 2);
        title('Result of A* Global Path Planner');

        % Apply B-Spline algorithm to smooth the generate path from A*
        smoothedPath = bsplinesmoothing(Map, path);
        plot(smoothedPath(:,2), smoothedPath(:,1), 'b--o', 'LineWidth', 2, 'MarkerSize', 2);
        title('Smoothed Path with B-spline Algorithm');
    else
        % If a path is not found, display 'Path not found' message in the command window
        disp(['Path from start pose (' num2str(start(1)) ', ' num2str(start(2)) ') to goal pose (' num2str(goal(1)) ', ' num2str(goal(2)) ') not found !']);
    end

    hold off;
end

%%
function path = astarpathplanner(Map, start, goal)
    % Initialize the open and closed lists as struct data type consist of
    % node, g cost, f cost as well as parent node
    open = struct('node', {}, 'g_cost', {}, 'f_cost', {}, 'parent', {});
    closed = struct('node', {}, 'g_cost', {}, 'f_cost', {}, 'parent', {});

    % Add the start node as a new element to the end of the open list
    open(end + 1) = struct('node', start, 'g_cost', 0, 'f_cost', calculate_heuristic(start, goal), 'parent', []);

    while ~isempty(open)
        % Find the node with the lowest f score in the open list and put it
        % out from the open list to become a current node
        [~, min_idx] = min(arrayfun(@(x) x.f_cost, open));
        current_node = open(min_idx);
        open(min_idx) = [];

        % If the goal node is reached, return the path
        if isequal(current_node.node, goal)
            path = build_path(current_node);
            return;
        end

        % Add the current node as a new element to the end of the closed list
        closed(end + 1) = current_node;

        % Expand the current node's neighbors (using 8 nearest node poses so it
        % can move diagonally)
        for i = -1:1
            for j = -1:1
                neighbor = current_node.node + [i, j];

                % Check if the neighbor is within the map bounds
                if all(neighbor >= [1, 1]) && all(neighbor <= size(Map))
                    % Check if the neighbor is occupied or unoccupied by
                    % obstacles
                    if Map(neighbor(1), neighbor(2)) == 0
                        % Calculate the g cost to reach the neighbor (using
                        % norm function so that diagonal movement will be
                        % 1.4 and other movement will be 1)
                        g_cost = current_node.g_cost + norm([i, j]);

                        % Calculate the heuristic cost to reach the goal from the neighbor
                        h_cost = calculate_heuristic(neighbor, goal);

                        % Calculate the total cost (f cost) to reach the neighbor
                        f_cost = g_cost + h_cost;

                        % Check if the neighbor is in the closed list
                        closed_idx = find(arrayfun(@(x) isequal(x.node, neighbor), closed));

                        if isempty(closed_idx)
                            % Check if the neighbor is in the open list
                            open_idx = find(arrayfun(@(x) isequal(x.node, neighbor), open));

                            if isempty(open_idx)
                                % Add the neighbor as a new element to the end of the open list
                                open(end + 1) = struct('node', neighbor, 'g_cost', g_cost, 'f_cost', f_cost, 'parent', current_node);
                            else
                                % Update the neighbor's information in the
                                % open list if the f cost is lower than the previous/old f cost value
                                if f_cost < open(open_idx).f_cost
                                    open(open_idx) = struct('node', neighbor, 'g_cost', g_cost, 'f_cost', f_cost, 'parent', current_node);
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    % Path not found will give an empty list
    path = [];
end

%%
% Function to build the path from the goal node to the start node
function path = build_path(node)
    path = [];
    
    % Loop to check the current node and its parent node
    while ~isempty(node.parent)
        path = [node.node; path];
        node = node.parent;
    end
    % Add the start node into the path list
    path = [node.node; path];
end

%%
% Function to calculate the heuristic cost to reach the goal node from the current
% node
function h_cost = calculate_heuristic(node, goal)
    % Calculate the Euclidean distance between the current node and the
    % goal node
    h_cost = norm(node - goal);
end

%%
function smoothedPath = bsplinesmoothing(Map, path)
    % B-spline interpolation to smooth the path
    t = 1:length(path);
    ts = 1:0.1:length(path);
    x = interp1(t, path(:, 2), ts, 'spline');
    y = interp1(t, path(:, 1), ts, 'spline');
    smoothedPath = [y' x'];

    % Check if the B-spline path is collision-free
    for i = 1:size(smoothedPath, 1)
        if Map(round(smoothedPath(i, 1)), round(smoothedPath(i, 2))) == 1
            disp('Path is not collision-free');
            return;
        end
    end
end