function draw_simulation(xt, zt, xh, zh, x, z, z_false, nt, nx, nz)
% DRAW_SIMULATION Plots true and estimated trajectories along with false alarms.
%
% Inputs:
%   xt       - ground truth state trajectories (cell array {k,t})
%   zt       - true measurements (cell array {k,t})
%   xh       - estimated state trajectories (cell array {k,t})
%   zh       - estimated measurements (cell array {k,t})
%   x        - 
%   z        -
%   z_false  - false alarms (cell array {t}{j})
%   T        - number of scans per target
%   nt       - number of targets
%   nx       - state dimension
%   nz       - measurement dimension

    T = size(xt,1);

    % Preallocate storage
    xv = zeros(T, nx, nt);
    zv = zeros(T, nz, nt);
    xhv = zeros(T, nx, nt);
    zhv = zeros(T, nz, nt);

    % Convert cell arrays to matrices
    for t = 1:nt
        for k = 1:T
            xv(k,1:nx,t) = xt{k,t};
            zv(k,1:nz,t) = zt{k,t};
            xhv(k,1:nx,t) = xh{k,t};
            zhv(k,1:nz,t) = zh{k,t};
            xnv(k,1:nx,t) = x{k,t};
            znv(k,1:nz,t) = z{k,t};
        end
    end

    % Plot
    figure();
    hold on;
    hnd = zeros(nt,1);
    lbl = cell(nt,1);

    for t = 1:nt
        color = rand(1,3);  % generate new color for each target
        plot(zv(:,1,t), zv(:,2,t), 'Color', color);
        hnd(t) = plot(zhv(:,1,t), zhv(:,2,t), 'o', 'MarkerFaceColor', color);
        noisy(t) = plot(znv(:,1,t), znv(:,2,t), 'x', 'MarkerFaceColor', color, 'MarkerSize', 10);
        lbl{t} = sprintf('Target %d', t);

        % Mark start and end
        plot(zv(1,1,t), zv(1,2,t), '^', 'Color', color, 'MarkerSize', 8, 'MarkerFaceColor', color);  % start
        plot(zv(end,1,t), zv(end,2,t), 'v', 'Color', color, 'MarkerSize', 8, 'MarkerFaceColor', color);  % end'

        % Mark last estimated measurement
        plot(zhv(end,1,t), zhv(end,2,t), 'kp', 'MarkerSize', 10, 'MarkerFaceColor', 'k');

    end

    % Draw false alarms on each frame
    for t = 1:length(z_false)
        for j = 1:length(z_false{t})
            if ~isempty(z_false{t}{j})
                plot(z_false{t}{j}(1), z_false{t}{j}(2), 'rx', 'LineWidth', 1.5);
            end
        end
    end

    % Final plot settings
    legend(hnd, lbl);
    set(gca, 'FontSize', 12);
    title('State-space', 'FontSize', 14);
    xlabel('Coordinate X (m)', 'FontSize', 14);
    ylabel('Coordinate Y (m)', 'FontSize', 14);
    grid on;

    legend(hnd,lbl);
    set(gca,'FontSize',12);
    title('State-space','FontSize',14);
    xlabel('Coordinate X (m)','FontSize',14);
    ylabel('Coordinate Y (m)','FontSize',14);

end