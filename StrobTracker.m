classdef StrobTracker < handle
    % STROBTRACKER  Plots strobe size & false marks from one queue,
    % and true measurements from another queue, for multiple trajectories.

    properties (Access = public)
        axesGrid   % numTraj×2 array of axes handles: col1 = X, col2 = Y
    end

    methods
        function obj = StrobTracker(measQueue, trueQueue, numTraj)
            % Constructor: create numTraj×2 subplots, legend, and subscribe
            %
            % measQueue : DataQueue sending strobe+center+false-mark events
            % trueQueue : DataQueue sending true-measurement events
            % numTraj   : number of trajectories (rows)

            figure('Name','Strobe Tracker','NumberTitle','off');
            obj.axesGrid = gobjects(numTraj,2);

            for idx = 1:numTraj
                % X-subplot
                axX = subplot(numTraj,2,(idx-1)*2+1); hold(axX,'on');
                xlabel(axX,'Measurement Index');
                ylabel(axX,'X Coordinate');
                title(axX,sprintf('Trajectory %d: X Strobe', idx));
                obj.axesGrid(idx,1) = axX;

                % Y-subplot
                axY = subplot(numTraj,2,(idx-1)*2+2); hold(axY,'on');
                xlabel(axY,'Measurement Index');
                ylabel(axY,'Y Coordinate');
                title(axY,sprintf('Trajectory %d: Y Strobe', idx));
                obj.axesGrid(idx,2) = axY;
            end

            % Shared legend on first X-subplot
            ax0 = obj.axesGrid(1,1);
            axes(ax0);  hold(ax0,'on');
            h1 = plot(NaN,NaN,'-','Color','green','LineWidth',2);
            h2 = plot(NaN,NaN,'x','Color','red','MarkerSize',8);
            h3 = plot(NaN,NaN,'.','Color','blue','MarkerSize',8);
            h4 = plot(NaN,NaN,'.','Color','red','MarkerSize',8);
            lg = legend(ax0, [h1,h2,h3,h4], ...
                   {'Strobe size','Center','True measurement','Marks'}, ...
                   'Location','best');
            lg.AutoUpdate = 'off';

            % Subscribe to measurement events (strobe/center/false)
            afterEach(measQueue, @obj.onMeasurement);

            % Subscribe to true-measurement events
            afterEach(trueQueue, @obj.onTrueMeasurement);
        end
    end

    methods (Access = private)
        function onMeasurement(obj, ev)
            % ev.trajectory – trajectory index
            % ev.step       – measurement index
            % ev.center     – [x; xdot; y; ydot]
            % ev.size       – [size_x; size_y]
            % ev.fm         – M×2 false marks [fx, fy]

            tr   = ev.trajectory;
            i    = ev.step;
            ctr  = ev.center;
            hs   = ev.size;
            fm   = ev.fm;

            % X subplot
            ax = obj.axesGrid(tr,1); axes(ax); hold(ax,'on');
            plot(ax, [i,i], [ctr(1)-hs(1), ctr(1)+hs(1)], ...
                 '-','Color','green','LineWidth',2);
            plot(ax, i, ctr(1), 'x','Color','red','MarkerSize',8);
            for k=1:size(fm,2)
                plot(ax, i, fm{k}(1), '.','Color','red','MarkerSize',8);
            end

            % Y subplot
            ax = obj.axesGrid(tr,2); axes(ax); hold(ax,'on');
            plot(ax, [i,i], [ctr(2)-hs(2), ctr(2)+hs(2)], ...
                 '-','Color','green','LineWidth',2);
            plot(ax, i, ctr(2),'x','Color','red','MarkerSize',8);
            for k=1:size(fm,2)
                plot(ax, i, fm{k}(2),'.','Color','red','MarkerSize',8);
            end
        end

        function onTrueMeasurement(obj, ev)
            % ev.trajectory – trajectory index
            % ev.step       – measurement index
            % ev.z          – 2×1 true measurement [zx; zy]

            tr  = ev.trajectory;
            i   = ev.step;
            zm  = ev.z;

            % X subplot: magenta star
            ax = obj.axesGrid(tr,1); axes(ax); hold(ax,'on');
            plot(ax, i, zm(1), '.','Color','magenta','MarkerSize',8);

            % Y subplot: magenta star
            ax = obj.axesGrid(tr,2); axes(ax); hold(ax,'on');
            plot(ax, i, zm(2), '.','Color','magenta','MarkerSize',8);
        end
    end
end