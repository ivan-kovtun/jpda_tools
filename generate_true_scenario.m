function [xt, cov_x] = generate_true_scenario(scenario_type, nt, T, Q)
% GENERATE_TRUE_SCENARIO
%   Generates true state trajectories based on the selected scenario type.
%
% Inputs:
%   scenario_type – string identifier of the scenario (e.g., 'narrow_tracks')
%   obj           – simulator object (for event notification if needed)
%   nt, T         – number of targets and time steps
%   Q             – process noise covariance
%   x0            – starting x-coordinate(s)
%   y_min, y_max  – y range
%
% Outputs:
%   xt            – T×nt cell of true state vectors
%   cov_x         – T×nt cell of process noise covariances (all = Q)

    

    switch lower(scenario_type)
        case 'narrow_tracks'

            y0 = 50;
            y_delta = 60;
            x0 = -90;

            [xt, cov_x] = generate_narrow_tracks_scenario(nt, T, Q, x0, y0, y_delta);
        otherwise
            error('Unsupported scenario type: %s', scenario_type);
    end
end