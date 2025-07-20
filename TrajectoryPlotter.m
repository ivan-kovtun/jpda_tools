classdef TrajectoryPlotter < handle
    properties (Access=private)
        ax            % axes handle
        colors        % map: trajectory → RGB color
        hTrue         % map: trajectory → line handle for true path
    end
    properties (Access=public)
        fig           % handle to the figure window
    end

    methods
        function obj = TrajectoryPlotter(publisher)
            % Create figure & axes
            obj.fig = figure('Name','Trajectory Plotter', ...
                             'NumberTitle','off', ...
                             'CloseRequestFcn', @(src,evt) delete(src), ...
                             'Visible','on', ...
                             'WindowStyle','normal');
            obj.ax = axes; 
            hold(obj.ax,'on');
            xlabel(obj.ax,'X'); 
            ylabel(obj.ax,'Y');
            title(obj.ax,'True & Measured Trajectories');

            % prepare maps
            obj.colors = containers.Map('KeyType','double','ValueType','any');
            obj.hTrue  = containers.Map('KeyType','double','ValueType','any');

            % Subscribe to the publisher’s NewMeasurement event
            addlistener(publisher, 'NewMeasurement', @obj.onNewMeasurement);
        end

        function onNewMeasurement(obj, ~, ev)
            traj = ev.trajectory;
            step = ev.step;
            tau  = ev.tau;
            x    = ev.x;   % [x; xdot; y; ydot]
            z    = ev.z;   % [z_x; z_y]

            % Assign a color & line handle on first sighting
            if ~isKey(obj.colors, traj)
                c = rand(1,3);
                obj.colors(traj) = c;
                h = plot(obj.ax, NaN, NaN, '-', ...
                         'Color',c, 'LineWidth',1.5, ...
                         'DisplayName',sprintf('Target %d',traj));
                obj.hTrue(traj) = h;
                legend(obj.ax,'show');
            else
                c = obj.colors(traj);
                h = obj.hTrue(traj);
            end

            % 1) Extend true‐trajectory line
            Xd = get(h,'XData'); 
            Yd = get(h,'YData');
            set(h, 'XData',[Xd, x(1)], 'YData',[Yd, x(3)]);

            % 2) Plot this noisy measurement
            plot(obj.ax, z(1), z(2), 'o', ...
                 'Color',c, 'MarkerFaceColor',c, ...
                    'HandleVisibility','off');

            % % 3) Label it with the measurement index
            % text(z(1), z(2), sprintf('%d',tau), ...
            %      'VerticalAlignment','bottom', ...
            %      'HorizontalAlignment','right', ...
            %      'Color',c, ...
            %      'FontSize',9, ...
            %      'HandleVisibility','off');

            % % Label for center
            % text(x(1)-0.8, x(2)-0.8, sprintf('%d', step), ...
            %      'Color', 'r', ...
            %      'FontSize', 9, ...
            %      'HandleVisibility','off');

            % 4) If tau==1, mark as start
            if tau == 1
                plot(obj.ax, z(1), z(2), '^', ...
                     'Color',c, 'MarkerSize',8, 'MarkerFaceColor',c, ...
                        'HandleVisibility','off');
            end

            legend(obj.ax, 'show');

            drawnow;
        end
    end
end