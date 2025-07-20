function [x, z, zt, z_false, cov_z, false_targets] = realize_scenario(...
    xt, cov_x, R, Rf, obs, gen_sys_noise, gen_obs_noise, lambda, box_size)

    [T, nt] = size(xt);
    nz = size(box_size, 1);

    x = cell(T, nt);
    z = cell(T, nt);
    zt = cell(T, nt);
    cov_z = cell(T, nt);
    z_false = cell(T, 1);
    false_targets = 0;

    for ti = 1:nt
        for k = 1:T
            % Add process noise
            w = gen_sys_noise(cov_x{k,ti});
            x{k,ti} = xt{k,ti} + w;

            % Select Rk
            if k < T
                Rk = R;
                vk = gen_obs_noise(R);
            else
                Rk = Rf;
                vk = gen_obs_noise(Rf);
            end

            % Get true measurement
            [z_true, S_true, ~] = obs(xt{k,ti}, cov_x{k,ti}, zeros(nz,1), zeros(nz));
            zt{k,ti} = z_true;
            cov_z{k,ti} = S_true;

            % Get noisy measurement
            [z_noisy, S_noisy, ~] = obs(xt{k,ti}, cov_x{k,ti}, vk, Rk);
            z{k,ti} = z_noisy;
            cov_z{k,ti} = S_noisy;
        end
    end

    % Clutter
    for k = 1:T
        nf = poissrnd(lambda);
        false_targets = false_targets + nf;
        z_false{k} = arrayfun(@(j) (rand(nz,1)-0.5).*box_size, 1:nf, 'UniformOutput', false);
    end
end