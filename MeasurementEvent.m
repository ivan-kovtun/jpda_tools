classdef MeasurementEvent < event.EventData
    properties
        trajectory   % integer ID of the target
        tau          % time‐step index
        x            % 4×1 true state [x; xdot; y; ydot]
        z            % 2×1 measured [z_x; z_y]
    end
    methods
        function obj = MeasurementEvent(traj, tau, x, z)
            obj.trajectory = traj;
            obj.tau        = tau;
            obj.x          = x;
            obj.z          = z;
        end
    end
end