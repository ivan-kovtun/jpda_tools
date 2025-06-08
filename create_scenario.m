function [sys_f, obs_f, x, xt, cov_x, cov_sys, z, z_false, zt, cov_z, cov_obs, false_targets] = create_scenario(nt, nx, nz, T, Q0, Q, R, Rf, sys, obs, ap, bp, av, bv, lambda, box_size, gen_sys_noise, gen_obs_noise)
    % createScenario - Simulates the system and generates the scenario for multiple targets.
    %
    % Inputs:
    %   nt            - Number of targets
    %   nx            - Number of states
    %   nz            - Number of observations
    %   T             - Number of time steps
    %   Q0            - Initial process noise covariance matrix
    %   Q             - Process noise covariance matrix
    %   R             - Measurement noise covariance matrix
    %   Rf            - Final measurement noise covariance matrix
    %   sys           - State transition function handle
    %   obs           - Observation function handle
    %   ap, bp        - Range for initial position
    %   av, bv        - Range for initial velocity
    %   lambda        - Spatial density of false measurements (clutter density)
    %   box_size      - Size of the area for false alarms
    %   gen_sys_noise - Function handle to generate process noise
    %   gen_obs_noise - Function handle to generate observation noise
    %
    % Outputs:
    %   sys_f         - Cell array of state transition functions for each target
    %   obs_f         - Cell array of observation functions for each target
    %   x             - Cell array of simulated states for each target
    %   xt            - Cell array of true states for each target
    %   cov_x         - Cell array of state covariance matrices
    %   cov_sys       - Cell array of process noise covariance matrices
    %   z             - Cell array of simulated observations for each target
    %   z_false       - Cell array of false alarms (clutter) for each time step
    %   zt            - Cell array of true observations for each target
    %   cov_z         - Cell array of observation covariance matrices
    %   cov_obs       - Cell array of measurement noise covariance matrices
    %   false_targets - Total number of false alarms generated
    
    % Initialize outputs
    sys_f = cell(1, nt);
    obs_f = cell(1, nt);
    x = cell(T, nt);
    xt = cell(T, nt);
    cov_x = cell(T, nt);
    cov_sys = cell(1, nt);
    z = cell(T, nt);
    z_false = cell(T, 1);
    zt = cell(T, nt);
    cov_z = cell(T, nt);
    cov_obs = cell(1, nt);
    false_targets = 0;
    
    % Simulate system for all targets
    for t = 1:nt
        % Assign the state and output functions to the target
        sys_f{t} = sys;
        obs_f{t} = obs;
    
        % Assign the process and measurement noise covariances to the target
        cov_sys{t} = Q;
        cov_obs{t} = R;
    
        % Initialize
        u = zeros(nx, T);
        v = zeros(nz, T);
        u(:, 1) = gen_sys_noise(Q0);     % initial process noise
        v(:, 1) = gen_obs_noise(R);      % initial observation noise
    
        x0 = zeros(nx, 1);
        x0([1 3], 1) = ap + (bp - ap) .* rand(2, 1); % initial position for target t
        x0([2 4], 1) = av + (bv - av) .* rand(2, 1); % initial velocity for target t
        x{1, t} = x0 + gen_sys_noise(Q);
        z{1, t} = obs(x0, Q0, v(:, 1), R);
    
        % True state and observation
        xt{1, t} = x{1, t};
        zt{1, t} = obs(x{1, t}, zeros(nx, nx), zeros(nz, 1), zeros(nz, nz));
    
        for k = 2:T
            u(:, k) = gen_sys_noise(Q);     % process noise
            if k < T
                v(:, k) = gen_obs_noise(R);
            else
                v(:, k) = gen_obs_noise(Rf);
            end
            x{k, t} = sys_f{t}(x{k - 1, t}, zeros(nx, nx), u(:, k), zeros(nx, nx));   % simulate state
            z{k, t} = obs_f{t}(x{k, t}, zeros(nx, nx), v(:, k), zeros(nz, nz));     % simulate observation
    
            % Add false alarms (clutter) to z_false{k}
            nf = poissrnd(lambda);  % number of false alarms
            false_targets = false_targets + nf;
            z_false{k} = cell(1, nf);
            for jj = 1:nf
                z_false{k}{jj} = (rand(nz, 1) - 0.5) .* box_size;  % false measurement
            end
    
            % True state and observation (without noise)
            xt{k, t} = x{k, t};
            zt{k, t} = obs_f{t}(x{k, t}, zeros(nx, nx), zeros(nz, 1), zeros(nz, nz));
        end
    end
    end