function [xt, cov_x] = generate_narrow_tracks_scenario(...
    nt, T, Q, x0, y0, y_delta)

    xt    = cell(T, nt);
    cov_x = cell(T, nt);

    t1 = 2; t2 = 1; t3 = 4; t4 = 1; t5 = 2;
    st = t1 + t2 + t3 + t4 + t5;
    T1 = round(T * t1 / st);
    T2 = round(T * t2 / st);
    T3 = round(T * t3 / st);
    T4 = round(T * t4 / st);
    T5 = T - T1 - T2 - T3 - T4;

    angle0 = [-20 * 2 * pi / 360, 20 * 2 * pi / 360];  
    v = 140 / (T - 1);
    x0 = x0 - 20 + 3 * rand(1, nt);
    y0 = y0 + [0, -y_delta];

    x_traj = zeros(nt, T);
    y_traj = zeros(nt, T);
    vx_traj = zeros(nt, T);
    vy_traj = zeros(nt, T);

    for ti = 1:nt
        angle = angle0(ti);
        x_traj(ti, 1) = x0(ti);
        y_traj(ti, 1) = y0(ti);

        % Segment 1
        for k = 2:T1
            x_traj(ti, k) = x_traj(ti, k-1) + v * cos(angle);
            y_traj(ti, k) = y_traj(ti, k-1) + v * sin(angle);
        end

        % Segment 2 (converge)
        for k = T1+1:T1+T2
            frac = (k-T1-1)/(T2-1);
            new_angle = angle * (1 - frac);
            x_traj(ti, k) = x_traj(ti, k-1) + v * cos(new_angle);
            y_traj(ti, k) = y_traj(ti, k-1) + v * sin(new_angle);
        end

        % Segment 3 (cruise)
        for k = T1+T2+1:T1+T2+T3
            x_traj(ti, k) = x_traj(ti, k-1) + v;
            y_traj(ti, k) = y_traj(ti, k-1);
        end

        % Segment 4 (diverge)
        for k = T1+T2+T3+1:T1+T2+T3+T4
            frac = (k-T1-T2-T3-1)/(T4-1);
            new_angle = -angle * frac;
            x_traj(ti, k) = x_traj(ti, k-1) + v * cos(new_angle);
            y_traj(ti, k) = y_traj(ti, k-1) + v * sin(new_angle);
        end

        % Segment 5
        for k = T1+T2+T3+T4+1:T
            x_traj(ti, k) = x_traj(ti, k-1) + v * cos(-angle);
            y_traj(ti, k) = y_traj(ti, k-1) + v * sin(-angle);
        end
    end

    % Fill true state trajectories and covariances
    for ti = 1:nt
        for k = 1:T
            if k > 1 
                vx = x_traj(ti, k) - x_traj(ti, k-1);
                vy = y_traj(ti, k) - y_traj(ti, k-1);
            else
                vx = 0;
                vy = 0;
            end
            
            xt{k, ti} = [x_traj(ti, k); vx; y_traj(ti, k); vy];
            cov_x{k, ti} = Q;
        end
    end
end