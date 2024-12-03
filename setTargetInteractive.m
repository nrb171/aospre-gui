classdef setTargetInteractive < handle
    properties 
        Targets = struct(...
            'map',[], ...    % assigned after region finding. 2D map. Values: 0=no target, integer = targetID.
            'weight',[], ... % default 10, can be set
            'dwellTime', [], ... % default 0 (no limit), can be set.
            'seedLocations', []);% same size as map. Best guess for target locations. Can be set via API or interactive app. 
        environment = []
        fig         % UIFigure handle for interactive target acquisition
        

    end

    properties (Access = private)
        ax          % UIAxes handle 
        pointEF     
        lineEF      
        polyEF      
        nIterationsEF
        nIterationsLabel
        targetID = 1; % starting value for targets (in case multiple are assigned)
        pointBtn
        lineBtn
        polyBtn
        doneBtn
        dwellTimeLabel
        dwellTime
        targetWeight
        targetWeightLabel

    end

    methods (Access = private)
        function obj = updatePlot(obj, backgroundBool, targetBool)
            if backgroundBool 
                hold(obj.ax, 'off')
                pcolor(obj.ax, obj.environment)
                shading(obj.ax, 'flat')
            end
            if targetBool 
                hold(obj.ax, 'on')
                uc = unique(obj.Targets.map);
                uc = uc(uc>0);
                cmap = hsv(numel(uc));
                for tID = 1:numel(uc)
                    contour(obj.ax, obj.Targets.map == tID, 1, 'LineColor', cmap(tID,:));
                end
            end
        end
        
        function obj = drawROI(obj, option)

            % container for all draw ROI calls. 
            % option = 'point': draw points until 'Esc' is pressed.
            % option = 'line': draw a polyline until 'Enter' is pressed.
            % option = 'poly': draw a closed polygon around target until
                % polygon is closed.

            switch option
                case "point"
                    % set mask and roi arrays before entry into loop
                    roi = {[]};
                    mask = zeros(([numel(obj.ax.Children(end).YData), numel(obj.ax.Children(end).XData)]));
                    % force entry into loop

                    
                    while 1
                        roi{end+1} = drawpoint(obj.ax);
                        try
                            mask(round(roi{end}.Position(2)), round(roi{end}.Position(1))) = 1;
                        catch ME
                            break
                        end
                    end
                    mask = imfilter(mask, fspecial("disk", 3));
                    cb = obj.pointEF.Value;
                case "line"
                    
                    roi = drawpolyline(obj.ax);
                    mask = createMask(roi, numel(obj.ax.Children(end).YData), numel(obj.ax.Children(end).XData));
                    cb = obj.lineEF.Value;
        
                case "poly"
                    roi = drawpolygon(obj.ax);
                    mask = createMask(roi, numel(obj.ax.Children(end).YData), numel(obj.ax.Children(end).XData));
                    cb = obj.polyEF.Value;
            end

            
            obj.Targets.seedLocations = mask;
            
            obj = obj.getTarget('ContractionBias', cb, 'MaxIterations', obj.nIterationsEF.Value);

            obj = obj.updatePlot(1,1);
            
            clear roi
            %uiwait(obj.fig)
        end % end drawROI()


        
    end

    
    methods (Access = public)
        function obj = getTarget(obj, varargin)
            % set done button to red until finalized.
            obj.doneBtn.BackgroundColor = [1 0.5 0.5];

            p = inputParser();
            addParameter(p, 'ContractionBias', 0)
            addParameter(p, 'MaxIterations', 200)
            parse(p, varargin{:});
            
            if isempty(obj.Targets.map)
                obj.Targets.map = zeros(size(obj.environment));
            end
            bw = activecontour(...
                obj.environment, ...
                obj.Targets.seedLocations, ...
                p.Results.MaxIterations, ...
                "Chan-vese", ...
                'ContractionBias', p.Results.ContractionBias);
            
            obj.Targets.map(obj.Targets.map == obj.targetID | bw) = obj.targetID;

        end % end getTarget()

        function obj = clearTarget(obj)
            obj.Targets.map = zeros(size(obj.Targets.map));
            obj = obj.updatePlot(1,0);
        end

        function obj = eraseTarget(obj)
            % delete contiguous regions
            roi = {};
            while 1
                mask = zeros(([numel(obj.ax.Children(end).YData), numel(obj.ax.Children(end).XData)]));
                roi{end+1} = drawpoint(obj.ax);
                try
                    mask(round(roi{end}.Position(2)), round(roi{end}.Position(1))) = 1;
                    bw = imreconstruct(logical(mask)*999999, obj.Targets.map);
                    obj.Targets.map(bw & obj.Targets.map == obj.targetID) = 0;
                    obj.updatePlot(1, 1)
                catch ME
                    break
                end
            end
        end % end of eraseTarget()
        
        function obj = finalizeTarget(obj)
            obj.Targets.weight(obj.targetID) = obj.targetWeight.Value;
            obj.Targets.dwellTime(obj.targetID) = obj.dwellTime.Value;
            obj.targetID = obj.targetID + 1;
            obj.doneBtn.BackgroundColor = [0.5 1 0.5];


        end

    end % end public methods

    %% methods for component creation
    methods (Access = private)
        function obj = createComponents(obj)
            % Draw the UI and assign logic functions to items, update
            % target map.
            if isempty(obj.environment)
                error('Set obj.environment before obj.drawGui() is called')
            end

            %% UIFigure and UIAxes
            obj.fig = uifigure('Position',[500 500 500 275]*2);        
            obj.ax = uiaxes(obj.fig, "Position",[150, 10, 340, 255]*2);

            obj.updatePlot(1,0)
            
            

            %% target buttons (clear, erase, finalize, etc.)
            clearBtn = uibutton(obj.fig, ...
                "ButtonPushedFcn", @(src,event) obj.clearTarget(), ...
                'Text', 'Clear Target', ...
                'Tooltip', "Clear current target and all ROIs on the figure.", ...
                'Position', [150, 200, 100, 25]);
            eraseBtn = uibutton(obj.fig, ...
                "ButtonPushedFcn", @(src,event) obj.eraseTarget(), ...
                'Text', 'Erase Tool', ...
                'Tooltip', "Use point ROI to delete contiguous regions of the current target.", ...
                'Position', [30, 250, 100, 25]);
            obj.doneBtn = uibutton(obj.fig, ...
                "ButtonPushedFcn", @(src,event) close(obj.fig), ...
                'Text', 'Done/Exit', ...
                'Tooltip', "Close the interactive target selection window.", ...
                'Position', [30, 200, 100, 25], ...
                'BackgroundColor', [1 0.5 0.5]);
            finalizeBtn = uibutton(obj.fig, ...
                "ButtonPushedFcn", @(src,event) obj.finalizeTarget(), ...
                'Text', 'Finalize Target', ...
                'Tooltip',"Finish selecting this target and/or move to the next target.", ...
                'Position', [150, 250, 100, 25]);
        
            %% region growth settings.
            % uilabel(obj.fig, "Text",{'Growth Tendency', '(0,1]=shrink; [-1,0) growth'}, ...
            %     "Position", [150, 220, 150, 50], "WordWrap","on");
            obj.pointEF = uieditfield(obj.fig, 'numeric', ...
                "Limits",[-1 1], ...
                "Position", [150, 400, 100, 25], ...
                "Value", -0.2);
            obj.lineEF = uieditfield(obj.fig,"numeric", ...
                "Limits",[-1 1], ...
                "Position", [150, 350, 100, 25], ...
                "Value", -0.1);
            obj.polyEF = uieditfield(obj.fig,"numeric", ...
                "Limits",[-1 1], ...
                "Position", [150, 300, 100, 25], ...
                "Value", 0.1);
            

            %% Extra target settings
            obj.dwellTimeLabel = uilabel(obj.fig, ...
                "Text", "Target Dwell Time: ", ...
                'Position',[30, 50, 100, 25]);
            obj.dwellTime = uieditfield(obj.fig,"numeric", ...
                "Limits",[0 inf], ...
                "Position", [150, 50, 100, 25], ...
                "Tooltip", "How long (in secondfs) do you want to observe target (useful in multiple-target scenarios).([0, inf]; 0=default. Note: 0=inf)", ...
                "Value", 0);

            obj.targetWeightLabel = uilabel(obj.fig, ...
                "Text", "Target weight for flight optimization.", ...
                'Position',[30, 100, 100, 25]);
            obj.targetWeight = uieditfield(obj.fig,"numeric", ...
                "Limits",[1 50], ...
                "Position", [150, 100, 100, 25], ...
                "Tooltip", "How heavily do you want to weight this target? ([1,50]; 10=default)", ...
                "Value", 10);

            obj.nIterationsLabel = uilabel(obj.fig, ...
                "Text", "Max Iterations: ", ...
                'Position',[30, 150, 100, 25]);
            obj.nIterationsEF = uieditfield(obj.fig,"numeric", ...
                "Limits",[1 1000], ...
                "Position", [150, 150, 100, 25], ...
                "Tooltip", "How many iterations do you want the region growing algorithm to perform? ([1,1000]; 200=default)", ...
                "Value", 200);

            %% ROI selection buttons
            obj.pointBtn = uibutton(obj.fig, ...
                "ButtonPushedFcn", @(src,event) obj.drawROI("point"), ...
                'Text', 'Draw Point ROI', ...
                'Position', [30, 400, 100, 25]);
            obj.lineBtn = uibutton(obj.fig, ...
                "ButtonPushedFcn", @(src,event) obj.drawROI("line"), ...
                'Text', 'Draw Line ROI', ...
                'Position', [30, 350, 100, 25]);
            obj.polyBtn = uibutton(obj.fig, ...
                "ButtonPushedFcn", @(src,event) obj.drawROI("poly"), ...
                'Text', 'Draw Polygon ROI', ...
                'Position', [30, 300, 100, 25]);

        end
    end
    %% GUI draw constructor
    methods (Access = public)
        function obj = setTargetInteractive(environment)
            obj.environment = environment;
            obj = obj.createComponents();

            uiwait(obj.fig)
            

        end
    end
end