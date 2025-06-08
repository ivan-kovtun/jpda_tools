classdef SimulatorCrossedTracks < handle
    events
        NewMeasurement
    end

    methods

        function [sys_f, obs_f, x, xt, cov_x, cov_sys, z, z_false, zt, cov_z, cov_obs, false_targets] = create_scenario(...
                obj, nt, nz, T, Q, R, Rf, sys, obs, r0, lambda, box_size, gen_sys_noise, gen_obs_noise)
            % CREATE_SCENARIO
            % Simulates nt targets on a piecewise converge–cruise–diverge trajectory,
            % and generates both noisy and noise‑free measurements via the supplied obs().
            %
            % **Note**: we do *not* call the constant‑velocity sys() inside the main loop—
            %       the states come directly from the precomputed x_traj/y_traj.
            %
            % Motion pattern (total ΔX = 50 + 70 + 50 = 170 m over T‑1 steps):
            %   1) First 50 m in X:  targets converge from random Y∈[y_min,y_max]
            %                       into a 5 m‑spaced formation.
            %   2) Next 70 m in X:   targets cruise together at constant spacing.
            %   3) Last 50 m in X:   targets diverge back to their original Y’s.
            %
            % Inputs:
            %   nt            – number of targets
            %   nx            – state dimension (here: 4 = [x;xdot;y;ydot])
            %   nz            – measurement dimension (here: 2 = [x;y])
            %   T             – number of time steps
            %   Q0            – initial process noise covariance (unused for state)
            %   Q             – per‑step process noise covariance
            %   R             – measurement noise covariance (for steps 1..T‑1)
            %   Rf            – measurement noise covariance at final step T
            %   sys           – state transition handle (ignored here)
            %   obs           – measurement handle:
            %                     [z_k, S_k, W_k] = obs(x_k, P_k, v_k, R_k)
            %   x0            – common starting X coordinate
            %   y_min, y_max  – range for random initial Y coordinates
            %   lambda        – Poisson rate of clutter per time step
            %   box_size      – nz×1 vector, extent of clutter box in each dimension
            %   gen_sys_noise – w = gen_sys_noise(Q)
            %   gen_obs_noise – v = gen_obs_noise(R)
            %
            % Outputs:
            %   sys_f         – 1×nt cell array of sys handles (passed through)
            %   obs_f         – 1×nt cell array of obs handles
            %   x             – T×nt cell of noisy state vectors
            %   xt            – T×nt cell of true (noise‑free) state vectors
            %   cov_x         – T×nt cell of process noise covariances (all = Q)
            %   cov_sys       – 1×nt cell of Q for each target
            %   z             – T×nt cell of noisy measurements
            %   zt            – T×nt cell of true (noise‑free) measurements
            %   cov_z         – T×nt cell of measurement covariances
            %   cov_obs       – 1×nt cell of R for each target
            %   z_false       – T×1 cell of clutter lists
            %   false_targets – total number of clutter points generated

            % 0) Preallocate all outputs
            sys_f        = repmat({sys},  1, nt);
            obs_f        = repmat({obs},  1, nt);
            cov_sys      = repmat({Q},    1, nt);
            cov_obs      = repmat({R},    1, nt);
            x            = cell(T, nt);
            xt           = cell(T, nt);
            cov_x        = cell(T, nt);
            z            = cell(T, nt);
            zt           = cell(T, nt);
            cov_z        = cell(T, nt);
            z_false      = cell(T, 1);
            false_targets = 0;

            %% 1) Compute the piecewise trajectory in X


            angle_min = pi;  
            angle_max = 3*pi/2;

            alpha = (angle_max - angle_min)/(nt - 1);

            alpha_traj = angle_min : alpha : angle_max;

            % total distance
            v = 150 / (T - 1);


            %% 2) Compute initial and formation Y positions
            % y0       = y_min + (y_max - y_min).*rand(nt,1);   % nt×1 random starts
            y0 = r0 * sin(alpha_traj);
            x0 = r0 * cos(alpha_traj);


            %% 4) Build Y and Y‑velocity trajectories
            y_traj = zeros(nt, T);
            x_traj = zeros(nt, T);
            vx_traj = zeros(nt, T);
            vy_traj = zeros(nt, T);
            for ti = 1:nt
                y_start  = y0(ti);
                x_start  = x0(ti);
                
                k = 1;
                y_traj(ti,k) = y_start;
                x_traj(ti,k) = x_start;
                
                % — 1
                for k = 2:T
                    y_traj(ti,k) = y_traj(ti,k - 1) + v * sin(-alpha_traj(ti));
                    x_traj(ti,k) = x_traj(ti,k - 1) - v * cos(-alpha_traj(ti));
                    vy_traj(ti,k) = v * sin(-alpha_traj(ti));
                    vx_traj(ti,k) = v * cos(-alpha_traj(ti));
                end


            end

            %% 5) Main simulation loop: states, measurements, clutter
            for ti = 1:nt
                for k = 1:T
                    % --- 5.1 True (noise‑free) state vector
                    state_true = [ x_traj(ti,k);
                        vx_traj(ti,k);
                        y_traj(ti,k);
                        vy_traj(ti,k) ];
                    % state_true = [ y_traj(ti,k);
                    %     vy_traj(ti,k);
                    %     x_traj(ti,k);
                    %     vx_traj(ti,k) ];
                    xt{k,ti}    = state_true;
                    cov_x{k,ti} = Q;

                    % --- 5.2 Noisy state (add process noise)
                    w          = gen_sys_noise(Q);
                    x{k,ti}    = state_true + w;

                    % --- 5.3 Select measurement noise for this step
                    if k < T
                        Rk = R;
                        vk = gen_obs_noise(R);
                    else
                        Rk = Rf;
                        vk = gen_obs_noise(Rf);
                    end

                    % --- 5.4 Noise‑free measurement via obs()
                    [z_true, S_true, ~] = obs(state_true, cov_x{k,ti}, zeros(nz,1), zeros(nz));
                    zt{k,ti}            = z_true;
                    cov_z{k,ti}         = S_true;

                    % --- 5.5 Noisy measurement via obs()
                    [z_noisy, S_noisy, ~] = obs(x{k,ti}, cov_x{k,ti}, vk, Rk);
                    z{k,ti}               = z_noisy;
                    cov_z{k,ti}           = S_noisy;

                    notify(obj, 'NewMeasurement', ...
                        MeasurementEvent(ti, k, x{k,ti}, z{k,ti}));

                    % --- 5.6 Clutter (false alarms)
                    nf = poissrnd(lambda);
                    false_targets = false_targets + nf;
                    z_false{k} = arrayfun(@(j) (rand(nz,1)-0.5).*box_size, 1:nf, 'UniformOutput', false);
                end
            end
        end
    end
end