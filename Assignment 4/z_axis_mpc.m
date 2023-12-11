function j = z_axis_mpc(K,dt,p_0,v_0,a_0,pt,vt,at)

    % Define the maximum/minimum constraint for velocity, acceleration, and jerk
    vmax = 6.0; % Maximum velocity constraint
    vmin = 1.0; % Minimum velocity constraint
    amax = 3.0; % Maximum acceleration constraint
    amin = 1.0; % Minimum acceleration constraint
    jerkmax = 2.0; % Maximum jerk constraint
    jerkmin = 2.0; % Minimum jerk constraint
    
    % Construct the prediction matrix for position, velocity, and acceleration
    [Tp, Tv, Ta, Bp, Bv, Ba] = getPredictionMatrix(K, dt, p_0, v_0, a_0);

    % Define the weight that will be used in the cost function
    w1 = 100.0; % Weight for tracking error
    w2 = 1.0; % Weight for velocity soft constraint
    w3 = 1.0; % Weight for acceleration soft constraint
    w4 = 1.0; % Weight for jerk (hard constraint)
    w5 = 1.0; % Weight for applying soft constraint

    % Construct the optimization problem
    H = blkdiag(w4 * eye(K) + w1 * (Tp' * Tp) + w2 * (Tv' * Tv) + w3 * (Ta' * Ta), w5 * eye(K));
    F = [w1 * (Bp - pt)' * Tp + w2 * (Bv - vt)' * Tv + w3 * (Ba - at)' * Ta, zeros(1,K)];
    
    % Construct soft constraints
    A = [Tv, zeros(K);-Tv, -eye(K); Ta, zeros(K);-Ta, zeros(K); zeros(size(Ta)), -eye(K)];
    b = [vmax * ones(20,1) - Bv; vmin * ones(20,1) + Bv; amax * ones(20,1) - Ba; amin * ones(K,1) + Ba; zeros(K,1)];

    % Define hard constraints on jerk as additional restriction
    lb = -jerkmin * ones(K,1);
    ub = jerkmax * ones(K, 1);

    % Define linear equality constraints
    Aeq = []; % not use equality constraint for this case
    Beq = []; % not use equality constraint for this case
    
    % Solve the optimization problem using quadratic programming solver
    J = quadprog(H, F, A, b, Aeq, Beq, lb, ub);
    
    % The output of the desired jerk
    j = J(1);

end