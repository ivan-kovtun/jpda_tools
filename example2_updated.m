%% clear memory, screen, and close all figures
clear, clc, close all;
addpath('export_fig');
addpath('matlab-tree');

global FIG
global STROB_TRACKER

% simulator = SimulatorCrossedTracks();
simulator = SimulatorNarrowTracks();
% simulator = SimulatorSimple Tracks();
FIG = TrajectoryPlotter(simulator).fig;



%% Constant velocity model parameters
T = 1;
F2 = [1 T; 0 1];
% G2 = [T^2/2; T];
% Q2 = G2*G2';
% A_q = 10;

A_q = 100;

association_threshold = 7;

Q2 = A_q*[T^3/3 T^2/2; T^2/2 T];

Fkm1 = blkdiag(F2,F2);
Hk = [1 0 0 0; 0 0 1 0];

%% Process equation for targets x[k] = sys_f(x[k-1], P[k-1], u[k-1], Q[k-1]);
nt = 02; % number of targets
nx = 04; % number of states
sys = @state_eq2;

show_strobes = false;

%% Prepare strobe tracker
measure_queue = [];
true_queue = [];
if show_strobes == true
    measure_queue = parallel.pool.DataQueue;
    true_queue = parallel.pool.DataQueue;
    tracker = StrobTracker(measure_queue, true_queue, nt);    % two trajectories
    STROB_TRACKER = tracker.axesGrid;
end

%% Observation equation for targets z[k] = obs_f(x[k], P[k], v[k], R[k]);
nz = 2; % number of observations
obs = @output_eq2;

%% PDF of process noise and noise generator function
nu = 4;
Q0 = 1e-2*blkdiag(Q2,Q2);
Q = 1e-3*blkdiag(Q2,Q2);
% Q0 = 1e-3*blkdiag(Q2,Q2);
% Q = 1e-4*blkdiag(Q2,Q2);
mux = zeros(nx,1);
p_sys_noise   = @(u, Qu) mvnpdf(u, mux, Qu);
% gen_sys_noise = @(Qu) mvnrnd(mux, Qu)';
gen_sys_noise = @(Qu) mvnrnd(mux, (Qu + Qu')/2)';

%% PDF of observation noise and noise generator function
nv = 2;
% R = 1e-3*eye(nv);
R = 1*eye(nv);
% Rf = 0.05*eye(nv);
Rf = R;
muz = zeros(nz,1);
p_obs_noise   = @(v, Rv) mvnpdf(v, muz, Rv);
gen_obs_noise = @(Rv) mvnrnd(muz, Rv)';

%% Initial state covariance matrix
% P0 = zeros(nx);
P0 = Q0;
S0 = Hk*(Fkm1*P0*Fkm1' + Q0)*Hk' + R;

%% Number of time steps
% T = 50;
% T = 35;
T = 75;

%% Separate memory space
sys_f = cell(1,nt);
obs_f = cell(1,nt);

x = cell(T,nt);
xt = cell(T,nt);
cov_x = cell(T,nt);
cov_sys = cell(1,nt);

z = cell(T,nt);
z_false = cell(T,1);
zt = cell(T,nt);
cov_z = cell(T,nt);
cov_obs = cell(1,nt);

u = zeros(nu,T);
v = zeros(nv,T);

ap = -5;
bp = +5;
av = -2;
bv = +2;

lambda = 1;
box_size = [220; 220];

y_min = 50;
y_max = -50;
x0 = -90;


filename = 'scenario_data_2.mat';
% filename = 'scenario_data_parallel_merged.mat';
% filename = 'scenario_data_parallel_merged_2.mat';
% filename = 'scenario_data_parallel_switched.mat';
filename = 'scenario_data_parallel_success.mat';

sys_f        = repmat({sys},  1, nt);
obs_f        = repmat({obs},  1, nt);
cov_sys      = repmat({Q},    1, nt);
cov_obs      = repmat({R},    1, nt);

% Нова функція для генерації істиних траекторій
[xt, cov_x] = generate_true_scenario('narrow_tracks', nt, T, Q);

% Лінійні траєкторії що перетинаються
% [sys_f, obs_f, x, xt, cov_x, cov_sys, z, z_false, zt, cov_z, cov_obs, false_targets] = simulator.create_scenario(nt, nz, T, Q, R, Rf, sys, obs, 75, lambda, box_size, gen_sys_noise, gen_obs_noise);

% Паралельні траєкторії
% [sys_f, obs_f, x, xt, cov_x, cov_sys, z, z_false, zt, cov_z, cov_obs, false_targets] = simulator.create_scenario(nt, nz, T, Q, R, Rf, sys, obs, x0, y_min, y_max, lambda, box_size, gen_sys_noise, gen_obs_noise);
                                                                                                                 
%[sys_f, obs_f, x, xt, cov_x, cov_sys, z, z_false, zt, cov_z, cov_obs, false_targets] = simulator.create_scenario(nt, nz, T, Q, R, Rf, sys, obs, x0, y_min, y_max, lambda, box_size, gen_sys_noise, gen_obs_noise);
%                                                                                                                 %nt, nz, T, Q, R, Rf, sys, obs, x0, y_min, y_max, lambda, box_size, gen_sys_noise, gen_obs_noise
% 
% save(filename, 'sys_f', 'obs_f', 'x', 'xt', 'cov_x', 'cov_sys', ...
%      'z', 'z_false', 'zt', 'cov_z', 'cov_obs', 'false_targets');
% 
% 
% load(filename, 'sys_f', 'obs_f', 'x', 'xt', 'cov_x', 'cov_sys', ...
%      'z', 'z_false', 'zt', 'cov_z', 'cov_obs', 'false_targets');

% [sys_f, obs_f, x, xt, cov_x, cov_sys, z, z_false, zt, cov_z, cov_obs, false_targets] = simulator.create_scenario(nt, nx, nz, T, Q0, Q, R, Rf, sys, obs, ap, bp, av, bv, lambda, box_size, gen_sys_noise, gen_obs_noise);

%% Allocate memory
xh = cell(T+1,nt);
% the size of zh = cell(T,nt);




% P = blkdiag([R(1,1), R(1,1)/T; R(1,1)/T, 2*R(1,1)/T^2], [R(2,2), R(2,2)/T; R(2,2)/T, 2*R(2,2)/T^2]);
% S = Hk*(Fkm1*P*Fkm1' + Q0)*Hk' + R;
P = P0;
S = S0;

%% Parameters
% Volume of validation region
gamma_ = chi2inv(0.99,nz);
lambda = 0.01;
epsilon = 1e-3;
% cnz = pi^(nz/2)/gamma(nz/2 + 1);
% Estimation of the total surveillance region
% (union of validation region for all targets)
% Vt = nt*cnz*sqrt(det(gamma_*(S)));
% lambda = 1/Vt;

params.k            = 1;                % initial iteration number
params.m            = nt;               % number of tracks
params.Nt           = nt;               % number of targets
params.cov_sys      = cov_sys;          % process noise covariance matrix
params.cov_obs      = cov_obs;          % measurement noise covariance matrix
params.PDt          = 0.95;             % detection probability of target
params.lambda       = lambda;           % spatial density of false measurements / clutter density
params.gamma        = chi2inv(0.99,nz); % gate threshold - probability (PG) - for confidence of 99% and nz degrees of freedom

% Nr Monte Carlo runs
Nr = 100;
NEES = zeros(T,1);
ERMS = zeros(T,1);

SE_x = zeros(T, nt);  % накопичення квадратів похибок по X
SE_y = zeros(T, nt);  % накопичення квадратів похибок по Y

Bias_x = zeros(T, nt);  % накопичення похибки по X
Bias_y = zeros(T, nt);  % накопичення похибки по Y

Var_x = zeros(T, nt); % аналітична дисперсія вимірювання 
Var_y = zeros(T, nt);


tic

total_metrics = struct('nCases', 0, ...
                          'nOK', 0, ...
                          'nSwitched', 0, ...
                          'nMerged', 0, ...
                          'nLost', 0, ...
                          'nResult', 0, ...
                          'CFT', 0)

for i = 1:Nr

    fprintf('Run = %d/%d\n',i,Nr);

    [x, z, zt, z_false, cov_z, false_targets] = realize_scenario(...
                xt, cov_x, R, Rf, obs, gen_sys_noise, gen_obs_noise, lambda, box_size);

    %% Set the initial state
    for t = 1:nt
        xh{1,t} = xt{1,t}; % як перший прогноз беремо другу істинну відмітку
        % zh{1,t} = obs_f{t}(x{1,t}, 0, 0, 0);
        zh{1, t} = zt{1,t};
        cov_x{1,t} = P0 + epsilon * eye(nx);
    
        % xh{1,t} = [z{2,t}(1,1); (z{2,t}(1,1)-z{1,t}(1,1))/T; z{2,t}(2,1); (z{2,t}(2,1)-z{1,t}(2,1))/T];
        % cov_x{1,t} = blkdiag([R(1,1), R(1,1)/T; R(1,1)/T, 2*R(1,1)/T^2], [R(2,2), R(2,2)/T; R(2,2)/T, 2*R(2,2)/T^2]);
    end


    % Estimate state
    for k = 2:T
        % fprintf('Iteration = %d/%d\n',k,T);

        % State estimation and filtered observation
        params.k = k;
        % z_all = [z(k-1,:), z_false{k-1}];
        z_all = [z(k,:), z_false{k}];
        [xh(k,:), cov_x(k,:), zh(k,:)] = jpda_filter(sys_f, obs_f, xh(k-1,:), cov_x(k-1,:), z_all, params, 'parametric', measure_queue);
        % [xh(k,:), cov_x(k,:), zh(k,:)] = jpda_filter(sys_f, obs_f, xh(k-1,:), cov_x(k-1,:), z_all, params, 'non-parametric', measure_queue);
        % [xh(k,:), cov_x(k,:), zh(k,:)] = jpda_filter(sys_f, obs_f, xh(k-1,:), cov_x(k-1,:), z_all, params, 'parametric', measure_queue);
        % [xh(k,:), cov_x(k,:), zh(k,:)] = jpda_filter(sys_f, obs_f, xh(k-1,:), cov_x(k-1,:), z_all, params, 'tree', measure_queue);
        % [xh(k,:), cov_x(k,:), zh(k,:)] = jpda_filter(sys_f, obs_f, xh(k-1,:), cov_x(k-1,:), z_all, params, 'lbp', measure_queue);

        % Computation of NEES
        NEESkt = 0;
        ERMSkt = 0;
        for t = 1:nt
            % NEESkt = NEESkt ...
            %     + 0.5*(xt{k,t} - xh{k,t})'*(cov_x{k,t}\(xt{k,t} - xh{k,t})) -0.5*nx*nt;
            NEESkt = NEESkt ...
                + (xt{k,t} - xh{k,t})'*(cov_x{k,t}\(xt{k,t} - xh{k,t}))/nx;
            ERMSkt = ERMSkt ...
                + sum((zt{k,t} - zh{k,t}).^2);


            dx = zt{k,t}(1) - zh{k,t}(1);   % похибка по X
            dy = zt{k,t}(2) - zh{k,t}(2);   % похибка по Y

            Bias_x(k,t) = Bias_x(k,t) + dx;
            Bias_y(k,t) = Bias_y(k,t) + dy;
            
            SE_x(k,t) = SE_x(k,t) + dx^2;
            SE_y(k,t) = SE_y(k,t) + dy^2;

            Pk = cov_x{k, t};    % 4x4
            Var_x(k, t) = Var_x(k, t) + Pk(1,1);  % дисперсія по X
            Var_y(k, t) = Var_y(k, t) + Pk(3,3);  % дисперсія по Y

        end
        NEES(k,1) = NEES(k,1) + NEESkt/nt;
        ERMS(k,1) = ERMS(k,1) + ERMSkt/nt;


    end

    % draw_simulation(xt, zt, xh, zh, x, z, z_false, nt, nx, nz);

    metrics = evaluate_jpda_metrics(xh, xt, T, nt, association_threshold, 2, T*0.95, T);

    total_metrics.nCases     = total_metrics.nCases     +   metrics.nCases;
    total_metrics.nOK        = total_metrics.nOK        +   metrics.nOK;
    total_metrics.nSwitched  = total_metrics.nSwitched  +   metrics.nSwitched;
    total_metrics.nMerged    = total_metrics.nMerged    +   metrics.nMerged;
    total_metrics.nLost      = total_metrics.nLost      +   metrics.nLost;
    total_metrics.nResult    = total_metrics.nResult    +   metrics.nResult;
    total_metrics.CFT        = total_metrics.CFT        +   metrics.CFT;

end

disp(total_metrics);

timerVal = tic;
toc

NEES = NEES/Nr;
ERMS = sqrt(ERMS/Nr);

MSE_x = SE_x / Nr;
MSE_y = SE_y / Nr;

RMSE_x = sqrt(MSE_x);
RMSE_y = sqrt(MSE_y);

Bias_x = Bias_x / Nr;
Bias_y = Bias_y / Nr;

Var_x = Var_x / Nr;
Var_y = Var_y / Nr;


figure;
titles = {'Target 1 - X', 'Target 1 - Y', 'Target 2 - X', 'Target 2 - Y'};

for t = 1:nt
    % --- X координата ---
    subplot(2,2,2*(t-1)+1);
    plot(1:T, Bias_x(:,t), '-',  'LineWidth', 1.5); hold on;
    plot(1:T, MSE_x(:,t),  '--', 'LineWidth', 1.5);
    plot(1:T, Var_x(:,t),  ':',  'LineWidth', 1.5);
    title(titles{2*(t-1)+1});
    xlabel('Time Step'); ylabel('X error stats');
    legend('Мат.сподівання ПВ', 'Квадр.ПВ', 'Аналіт.оцінка похибки');
    grid on;

    % --- Y координата ---
    subplot(2,2,2*(t-1)+2);
    plot(1:T, Bias_y(:,t), '-',  'LineWidth', 1.5); hold on;
    plot(1:T, MSE_y(:,t),  '--', 'LineWidth', 1.5);
    plot(1:T, Var_y(:,t),  ':',  'LineWidth', 1.5);
    title(titles{2*(t-1)+2});
    xlabel('Time Step'); ylabel('Y error stats');
    legend('Мат.сподівання ПВ', 'Квадр.ПВ', 'Аналіт.оцінка похибки');
    grid on;
end




% Ступені свободи (для nt цілей)
dof = nx * nt;
alpha = 0.05;

% Границі довірчого інтервалу (95%)
lower_bound = chi2inv(alpha/2, dof) / nt;
upper_bound = chi2inv(1 - alpha/2, dof) / nt;

% Вивести межі
disp(['Lower bound (95%): ', num2str(lower_bound)])
disp(['Upper bound (95%): ', num2str(upper_bound)])

% xv = zeros(T,nx,nt);
% zv = zeros(T,nz,nt);
% xhv = zeros(T,nx,nt);
% zhv = zeros(T,nz,nt);
% 
% for t = 1:nt
%     for k = 1:T
%         xv(k,1:nx,t) =  xt{k,t};
%         zv(k,1:nz,t) =  zt{k,t};
%         xhv(k,1:nx,t) =  xh{k,t};
%         zhv(k,1:nz,t) =  zh{k,t};
%     end
% end
% 
% % Plot
% figure()
% hnd = zeros(nt,1);
% color = rand(3,1);
% plot(zv(:,1,1), zv(:,2,1), 'Color', color); hold on;
% hnd(1) = plot(zhv(:,1,1), zhv(:,2,1), 'o', 'MarkerFaceColor', color);
% lbl = cell(nt,1);
% lbl{1,1} = 'Target 1';
% 
% % Mark start and end
% t = 1;
% plot(zv(1,1,t), zv(1,2,t), '^', 'Color', color, 'MarkerSize', 8, 'MarkerFaceColor', color);  % старт
% plot(zv(end,1,t), zv(end,2,t), 'v', 'Color', color, 'MarkerSize', 8, 'MarkerFaceColor', color);  % фініш
% 
% for t = 2:nt
%     color = rand(3,1);
%     plot(zv(:,1,t), zv(:,2,t), 'Color', color);
%     hnd(t) = plot(zhv(:,1,t), zhv(:,2,t), 'o', 'MarkerFaceColor', color);
%     lbl{t,1} = sprintf('Target %d', t);
% 
%     % Mark start and end
%     plot(zv(1,1,t), zv(1,2,t), '^', 'Color', color, 'MarkerSize', 8, 'MarkerFaceColor', color);  % старт
%     plot(zv(end,1,t), zv(end,2,t), 'v', 'Color', color, 'MarkerSize', 8, 'MarkerFaceColor', color);  % фініш
% 
% end
% 
% % --- Draw false alarms on frame 1 ---
% for t = 1:length(z_false)
%     for j = 1:length(z_false{t})
%         if ~isempty(z_false{t}{j})
%             plot(z_false{t}{j}(1), z_false{t}{j}(2), 'rx', 'LineWidth', 1.5);
%         end
%     end
% end



figure()
plot((1:T)',NEES);
ylabel('NEES','FontSize',14);
xlabel('Epoch','FontSize',14);
set(gca,'FontSize',12);
title(sprintf('Normalised Estimation Error Squared (NEES) - %d targets',nt),'FontSize',14);
grid on;

figure()
plot((1:T)',ERMS);
ylabel('RMSE (m)','FontSize',14);
xlabel('Epoch','FontSize',14);
set(gca,'FontSize',12);
title(sprintf('Root-Mean-Square Error (RMSE) - %d targets',nt),'FontSize',14);
grid on;

return;
