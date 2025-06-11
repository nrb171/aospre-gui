classdef aospreguiapp < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                        matlab.ui.Figure
        TabGroup                        matlab.ui.container.TabGroup
        SimulationDetailsTab            matlab.ui.container.Tab
        SimResolutionLabel              matlab.ui.control.Label
        SimulationResolutionEditField   matlab.ui.control.NumericEditField
        FileBrowserButton               matlab.ui.control.Button
        TimeGlobPattern                 matlab.ui.control.EditField
        FileTimePatternEditFieldLabel   matlab.ui.control.Label
        MessagesUI1                     matlab.ui.control.Label
        SimulationGlobPattern           matlab.ui.control.EditField
        FileSearchPatternEditFieldLabel  matlab.ui.control.Label
        FlightPlanningTab               matlab.ui.container.Tab
        MessagesUI2                     matlab.ui.control.Label
        UITable                         matlab.ui.control.Table
        AutomaticFlightTargetPlanningPanel  matlab.ui.container.Panel
        StartOrContinueFlight           matlab.ui.control.DropDown
        InteractiveFlightPlannerButton  matlab.ui.control.Button
        InteractiveTargetPlannerButton  matlab.ui.control.Button
        ManualFlightPlanningPanel       matlab.ui.container.Panel
        ManuallyaddflightlegButton      matlab.ui.control.Button
        EditflightlegverticesButton     matlab.ui.control.Button
        CalculateBackgroundMotionPanel  matlab.ui.container.Panel
        CalculateBackgroundMotionButton  matlab.ui.control.Button
        AOSPREDetailsTab                matlab.ui.container.Tab
        GridLayout4                     matlab.ui.container.GridLayout
        MessagesUI3                     matlab.ui.control.Label
        InitAOSPREButton                matlab.ui.control.Button
        SameasplottedtimeCheckBox       matlab.ui.control.CheckBox
        StartTimeSlider                 matlab.ui.control.Slider
        AOSPREPath                      matlab.ui.control.EditField
        pathtoAOSPRELabel               matlab.ui.control.Label
        AircraftPanelScanSettings       matlab.ui.container.Panel
        GridLayout                      matlab.ui.container.GridLayout
        AircraftVelocity                matlab.ui.control.NumericEditField
        AircraftVelmsEditFieldLabel     matlab.ui.control.Label
        AftEditField                    matlab.ui.control.NumericEditField
        AftEditFieldLabel               matlab.ui.control.Label
        ForeEditField                   matlab.ui.control.NumericEditField
        ForeEditFieldLabel              matlab.ui.control.Label
        MaximumTiltRotEditField         matlab.ui.control.NumericEditField
        MaximumTiltRotLabel             matlab.ui.control.Label
        AngularResolutionEditField      matlab.ui.control.NumericEditField
        AngularResolutionLabel          matlab.ui.control.Label
        Tree                            matlab.ui.container.CheckBoxTree
        LHSPortPanelNode                matlab.ui.container.TreeNode
        LHSRHI                          matlab.ui.container.TreeNode
        LHSPPI                          matlab.ui.container.TreeNode
        LHSFore                         matlab.ui.container.TreeNode
        AftNode_2                       matlab.ui.container.TreeNode
        RHSStarboardPanelNode           matlab.ui.container.TreeNode
        RHSRHI                          matlab.ui.container.TreeNode
        RHSPPI                          matlab.ui.container.TreeNode
        RHSFore                         matlab.ui.container.TreeNode
        AftNode                         matlab.ui.container.TreeNode
        TopPanelNode                    matlab.ui.container.TreeNode
        TOPPPI                          matlab.ui.container.TreeNode
        TOPFore                         matlab.ui.container.TreeNode
        BottomPanelNode                 matlab.ui.container.TreeNode
        BOTFore                         matlab.ui.container.TreeNode
        AftNode_3                       matlab.ui.container.TreeNode
        StartTimeEditField              matlab.ui.control.NumericEditField
        StartTimeEditFieldLabel         matlab.ui.control.Label
        AdvancedOptionsTab              matlab.ui.container.Tab
        ScanningOptionsPanel            matlab.ui.container.Panel
        ScanningOptionsGrid             matlab.ui.container.GridLayout
        skip_seconds_between_scans      matlab.ui.control.NumericEditField
        skip_seconds_between_scansLabel  matlab.ui.control.Label
        max_range_in_meters             matlab.ui.control.NumericEditField
        max_range_in_metersEditFieldLabel  matlab.ui.control.Label
        meters_to_center_of_first_gate  matlab.ui.control.NumericEditField
        meters_to_center_of_first_gateEditFieldLabel  matlab.ui.control.Label
        meters_between_gates            matlab.ui.control.EditField
        meters_between_gatesEditFieldLabel  matlab.ui.control.Label
        CRSIM_ConfigEditFieldLabel      matlab.ui.control.Label
        CRSIM_ConfigPath                matlab.ui.control.EditField
        beams_per_acquisition_time      matlab.ui.control.NumericEditField
        beams_per_acquisition_timeEditFieldLabel  matlab.ui.control.Label
        revisits_per_acquisition_time   matlab.ui.control.NumericEditField
        revisits_per_acquisition_timeEditFieldLabel  matlab.ui.control.Label
        pulses_per_pulse_set            matlab.ui.control.NumericEditField
        pulses_per_pulse_setEditFieldLabel  matlab.ui.control.Label
        pulse_repetition_frequency      matlab.ui.control.NumericEditField
        pulse_repetition_frequencyEditFieldLabel  matlab.ui.control.Label
        flight_level_waypoints_vert     matlab.ui.control.EditField
        flight_level_waypoints_vertEditFieldLabel  matlab.ui.control.Label
        flight_level_coordinate         matlab.ui.control.EditField
        flight_level_coordinateEditFieldLabel  matlab.ui.control.Label
        output_filename_format_string   matlab.ui.control.EditField
        outputformatstringLabel         matlab.ui.control.Label
        PlotDetailsPanel                matlab.ui.container.Panel
        GridLayout3                     matlab.ui.container.GridLayout
        WarpTargetsCheckBox             matlab.ui.control.CheckBox
        WarpFlightPathCheckBox          matlab.ui.control.CheckBox
        AltToPlotUnits                  matlab.ui.control.DropDown
        SimLevelLabel_2                 matlab.ui.control.Label
        AltToPlotText                   matlab.ui.control.NumericEditField
        VariableToPlot                  matlab.ui.control.DropDown
        PlottedVariableDropDownLabel    matlab.ui.control.Label
        SimLevelLabel                   matlab.ui.control.Label
        LevelToPlot                     matlab.ui.control.NumericEditField
        LevelToPlotSlider               matlab.ui.control.Slider
        TimeIndexSlider                 matlab.ui.control.Slider
        TimeIndexLabel                  matlab.ui.control.Label
        TimeIndexEditField              matlab.ui.control.NumericEditField
        UIAxes                          matlab.ui.control.UIAxes
        ContextMenu                     matlab.ui.container.ContextMenu
        ClearMenu                       matlab.ui.container.Menu
        ClearSelection                  matlab.ui.container.Menu
        LegRowMenu                      matlab.ui.container.Menu
        RefreshMenu                     matlab.ui.container.Menu
        RefreshSelectionMenu            matlab.ui.container.Menu
        AddnewLegMenu                   matlab.ui.container.Menu
        ContextMenu2                    matlab.ui.container.ContextMenu
        RefreshMenu_2                   matlab.ui.container.Menu
    end

    
    properties (Access = public)
        FileList        % list of filenames from dir() call.
        FileTimes       % times (in MATLAB durations) of the FileList files.
        PanelLocations  % structure of Rotation/tilt for each panel.
        ScanTables      % structure containing all of the scan tables.
        
        environment     % array of environment vertical slice.
        messages = {}        % cell array of messages to communicate app actions to user.
        verticalLevels
        var
        loadedVariable = ''
        xlat
        xlon
        deformationTransforms = {} % cell array containing deformation transforms between timesteps of FileList.
        deformationOriginTimeStep = 1 % time step at which to begin the deformations for flight and targets. 
        deformationLoadedBool = 0
    end
    properties (Access = public, Dependent)

    end
    properties (Access = public)
        flightPath          % roi polyline of flight path.     
        Targets = struct(...
            'map',[], ...    % assigned after region finding. 2D map. Values: 0=no target, integer = targetID.
            'weight',[], ... % default 10, can be set
            'dwellTime', [], ... % default 0 (no limit), can be set.
            'loadedBool', 0)    % bool to let the shutdown/startup functions to save/load the targets.

    end
    
    methods (Access = private)
        
        function [] = UpdateFigure(app)
            % function to redraw plan view
            

            
            %% load data if a different variable or timestep is selected.
            if ~strcmp([app.VariableToPlot.Value,app.FileList(app.TimeIndexEditField.Value).name], app.loadedVariable)
                app.updateMessages('Loading data... ')
                app.var = app.loadVar(app.VariableToPlot.Value);
                app.loadedVariable = [app.VariableToPlot.Value,app.FileList(app.TimeIndexEditField.Value).name];
            end

            %calculate the deformation of targets and flights at new
            %timestep.
            if app.deformationOriginTimeStep ~= app.TimeIndexEditField.Value
                
                
                if app.WarpTargetsCheckBox.Value | app.WarpFlightPathCheckBox.Value
                    timeStepLimits = [app.deformationOriginTimeStep, app.TimeIndexEditField.Value];
                    timeStepsBetween = min(timeStepLimits):max(timeStepLimits)-1;

                    
                    %{
                    D1 = zeros(size(app.deformationTransforms{1}, [1,2]));
                    D2 = zeros(size(app.deformationTransforms{1}, [1,2]));
                    
                    for iTime = timeStepsBetween
                        D1 = D1 - app.deformationTransforms{iTime}(:,:,1)*sign(timeStepLimits(2)-timeStepLimits(1));
                        D2 = D2 - app.deformationTransforms{iTime}(:,:,2)*sign(timeStepLimits(2)-timeStepLimits(1));
                    end

                    tform = cat(3,D1,D2);
                    %}

                    
                   
                    A = diag([1,1,1]);
                    offDiagInds = find(A~=1);
                    for iTime = timeStepsBetween
                        A(offDiagInds) = A(offDiagInds) + app.deformationTransforms{iTime}.A(offDiagInds);
                    end
                    
                    
                    %flip transform if the direction is backwards in time.
                    A(offDiagInds) = A(offDiagInds).*sign(timeStepLimits(2)-timeStepLimits(1));
                    
                    Atemp = A;
                    A(1,3) = Atemp(2,3);
                    A(2,3) = Atemp(1,3);
                    A(2,1) = Atemp(1,2);
                    A(1,2) = Atemp(2,1);
                    tform = affinetform2d(A);
                    
                end
                    

                if app.WarpTargetsCheckBox.Value
                    
                    warpedTargets = imwarp(app.Targets.map, tform, 'nearest', 'OutputView',imref2d(size(app.Targets.map)));
                    
                    %warpedTargets = imwarp(app.Targets.map, -round(permute(tform, [2,1,3])), 'nearest');
                    app.Targets.map = warpedTargets;
                end
                
                %keyboard
                if  app.WarpFlightPathCheckBox.Value
                    % transform table coordinates
                    for iRow = 1:size(app.UITable.Data, 1)
                        u1 = app.UITable.Data(iRow,1);
                        v1 = app.UITable.Data(iRow,2);
                        [x1,y1] = transformPointsForward(tform, u1, v1);
                        %x1=u1+tform(round(u1), round(v1),1);
                        %y1=v1+tform(round(u1), round(v1),2);

                        app.UITable.Data(iRow,1) = x1;
                        app.UITable.Data(iRow,2) = y1;
    
    
                        u2 = app.UITable.Data(iRow,3);
                        v2 = app.UITable.Data(iRow,4);
                        [x2,y2] = transformPointsForward(tform, u2, v2);
                        %x2=u2+tform(round(u2), round(v2),1);
                        %y2=v2+tform(round(u2), round(v2),2);
                        app.UITable.Data(iRow,3) = x2;
                        app.UITable.Data(iRow,4) = y2;
                    end

                end
                app.deformationOriginTimeStep = app.TimeIndexEditField.Value;
                
                

            end


            
            
            if isempty(app.xlat)
                app.updateMessages('Loading XLAT... ', 'append')
                app.xlat = app.loadVar('XLAT');
               
                ylim(app.UIAxes, [min(app.xlat(:))-0.5, max(app.xlat(:))+0.5])
            
                app.updateMessages('Loading XLONG... ', 'append')
                app.xlon = app.loadVar('XLONG');
                xlim(app.UIAxes, [min(app.xlon(:))-0.5, max(app.xlon(:))+0.5])
            
                if all(app.xlat(:) == 0) || all(app.xlon(:) == 0)
                    [xlonTemp, xlatTemp] = meshgrid(1:size(app.xlat,1), 1:size(app.xlat,2));
                    app.xlat = xlatTemp;
                    app.xlon = xlonTemp;
                end
            end

            app.updateMessages('Redrawing figure... ')
            


            
            varSlice = app.var(:,:,app.LevelToPlot.Value);
            app.environment = varSlice';

            %% draw map
            hold(app.UIAxes, 'off')

            nuimagesc(app.UIAxes, app.xlat, app.xlon, varSlice);
            hold(app.UIAxes, 'on')

            plotEarth(app.UIAxes)

            colorbar(app.UIAxes)
            clim(app.UIAxes, [prctile(varSlice(:), 1), prctile(varSlice(:), 99)+1e-3])
            
            title(app.UIAxes, [app.VariableToPlot.Value, ' (~', num2str(app.verticalLevels.(app.AltToPlotUnits.Value)(app.LevelToPlot.Value), '%.0f '), ' ', app.AltToPlotUnits.Value, ')'])

            %% draw path if one exists
            if ~isempty(app.UITable.Data)
                hold(app.UIAxes, 'on')

                [x,y] = app.getWarpedFlight(0);  

                lons = [];
                lats = [];
                for i = 1:numel(x)-1
                    

                    [lat1,lon1] = app.rowcolTolatlon(x(i), y(i));
                    [lat2,lon2] = app.rowcolTolatlon(x(i+1), y(i+1));

                    dlat = lat2-lat1;
                    dlon = lon2-lon1;

                    
                    scatter(app.UIAxes, lon1, lat1, 'r', 'filled')
                    scatter(app.UIAxes, lon2, lat2, 'r', 'filled')

                  
                    h=quiver(app.UIAxes,lon1,lat1,dlon,dlat,0, 'filled', 'r', 'LineWidth',2);
                    
                    lons(end+1) = lon1;
                    lats(end+1) = lat1;

                    
                    %set(h,'MaxHeadSize',1,'AutoScaleFactor',1);
                    
                end

                lons(end+1) = lon2;
                lats(end+1) = lat2;
                P = polybuffer([lons(:), lats(:)], 'lines', app.max_range_in_meters.Value/111000);
                p = plot(app.UIAxes, P);
                p.FaceAlpha = 0.1;
                p.FaceColor = 'r';


                %draw the time-warped path adding background motion
                if app.WarpFlightPathCheckBox.Value == 1
                    [x,y] = app.getWarpedFlight(1);
                    for i = 1:numel(x)-1
    
                        [lat1,lon1] = app.rowcolTolatlon(x(i), y(i));
                        [lat2,lon2] = app.rowcolTolatlon(x(i+1), y(i+1));
    
                        dlat = lat2-lat1;
                        dlon = lon2-lon1;
    
                        
                        scatter(app.UIAxes, lon1, lat1, 'm', 'filled')
                        scatter(app.UIAxes, lon2, lat2, 'm', 'filled')
                        h=quiver(app.UIAxes,lon1,lat1,dlon,dlat,0, 'filled', 'm', 'LineWidth',2);
                        
                        
                        %set(h,'MaxHeadSize',1,'AutoScaleFactor',1);
                        
                    end
                end
                
                

                
                %plot(app.UIAxes, x, y, '-o', 'Color', 'r', 'LineWidth', 1.5)
            end

            %% draw Targets.map (if available)
            if ~isempty(app.Targets.map)
                uTargets = unique(app.Targets.map);
                uTargets = uTargets(uTargets>0);
                hold(app.UIAxes, 'on')
                cmap = [...
                    251,154,153,
                    227,26,28,
                    253,191,111,
                    255,127,0,
                    202,178,214,
                    106,61,154,
                    255,255,153,
                    177,89,40]/255;
                for i = numel(uTargets)
                    contour(app.UIAxes, app.xlon', app.xlat',app.Targets.map == uTargets(i), 1, 'LineColor', cmap(i,:), 'LineWidth', 2)
                end
            end

            %% clear datatips
            
            graphicsObjects = {app.UIAxes.Children.Children};
            for i = 1:numel(graphicsObjects)
                if isa(graphicsObjects{i},'matlab.graphics.datatip.DataTip')
                    delete(graphicsObjects{i})
                end
            end
        end
        
        function [] = SetFigureControlsVisibility(app)
            app.PlottedVariableDropDownLabel.Visible = "on";
            app.LevelToPlotSlider.Visible = "on";
            app.SimLevelLabel.Visible = "on";
            app.LevelToPlot.Visible = "on";
            app.TimeIndexEditField.Visible = "on";
            app.TimeIndexLabel.Visible = "on";
            app.TimeIndexSlider.Visible = "on";
            app.VariableToPlot.Visible = "on";
            app.SimResolutionLabel.Visible = "on";
            app.SimulationResolutionEditField.Visible = "on";
            app.AltToPlotText.Visible = "on";
            app.AltToPlotUnits.Visible = "on";
            
        end
        
        function [] = UpdateSlidersAndRanges(app)
            %% set simulation level 
            % get temporary variable to set number of levels/times
            
            
            % only update the vertical levels/sliders the first time data
            % is loaded.

            if isinf(app.LevelToPlot.Limits(1)) && isinf(app.LevelToPlot.Limits(2))
                
                tempVar = ncread(...
                    [...    %filename
                        app.FileList(1).folder,'/', ...
                        app.FileList(1).name...
                    ], ...  %variable name
                    app.VariableToPlot.Value ...
                );

                app.updateMessages('Calculating time and altitude ranges... ')

                % get levels/times
                numVerticalLevels = size(tempVar,3);
                timeIndices = numel(app.FileList);
                
                % set ranges and temp values for sliders/edit fields
                app.TimeIndexEditField.Value = 1;
                app.TimeIndexSlider.Limits = [1, timeIndices];
                app.TimeIndexSlider.MinorTicks = ceil(linspace(1, timeIndices, min(numel(app.FileList),50)));
                app.TimeIndexSlider.MajorTicks = intersect(round(linspace(1, timeIndices, ceil(min(sqrt(timeIndices), 10)))), app.TimeIndexSlider.MinorTicks);
                app.TimeIndexSlider.Value = 1;
                
    
                [~,bestLevelGuess] = max(squeeze(std(std(tempVar))));
                
                app.LevelToPlotSlider.Value = bestLevelGuess;
                app.LevelToPlotSlider.Limits = [1, numVerticalLevels];
                app.LevelToPlotSlider.MinorTicks = 1:numVerticalLevels;
                app.LevelToPlotSlider.MajorTicks = intersect(ceil(linspace(1, numVerticalLevels, sqrt(numVerticalLevels))), app.LevelToPlotSlider.MinorTicks);

                

                

                % calculate the vertical levels
                % km
                PH = app.loadVar('PH');
                PHB = app.loadVar('PHB');
                app.verticalLevels.m = squeeze(nanmean((PH+PHB)/9.81, [1,2]));

                clear PHB PH

                % hPa
                P = app.loadVar('P');
                PB = app.loadVar('PB');
                app.verticalLevels.hPa = squeeze(nanmean((P+PB)/100, [1,2]));

                app.setAltitudeSlidersAndValues(bestLevelGuess)

                app.updateMessages('\color{green}done!', 'append')

            end
        
        end
        
        function [] = UpdateTableLogic(app)
            % try to update the tables based on some simple logic
            app.updateMessages('Recalculating table variables.')

            %% extract position from ROI polyline
            % add empty row if none are there
            if isempty(app.UITable.Data)
                app.UITable.Data(1,:) = [NaN, NaN, NaN, NaN, NaN, NaN];
            end
            
            
            %loop through roi and add to figure
            if nnz(isnan(app.UITable.Data)) == numel(app.UITable.Data)
                startIdx = 1;
            else
                startIdx = 2;
            end
            try
                for iData = startIdx:size(app.flightPath.Position,1)
                    if nnz(isnan(app.UITable.Data(end, 1:2)))==2
                        %first entry initializes start of leg
                        app.UITable.Data(end,1:2) = app.flightPath.Position(iData,:);
    
                    elseif nnz(isnan(app.UITable.Data(end, 3:4)))==2
                        % if a first entry already exists, add to end of
                        % leg
                        app.UITable.Data(end,3:4) = app.flightPath.Position(iData,:);
    
                    else
                        % create new leg and repeat above processes
                        app.UITable.Data(end+1,:) = [NaN, NaN, NaN, NaN, NaN, NaN];
                        app.UITable.Data(end,1:2) = app.UITable.Data(end-1,3:4);
                        app.UITable.Data(end,3:4) = app.flightPath.Position(iData,:);
                    end
                end
                
                app.flightPath.delete()
            catch ME

            end
                
            

            



            TAB = app.UITable.Data;
            for iRow = 1:size(app.UITable.Data,1)
                try
                    % set start of next leg to end of last leg
                    if iRow >1 && nnz(isnan(TAB(iRow, 1:2)))==2
                        TAB(iRow, 1) = TAB(iRow-1, 3);
                        TAB(iRow, 2) = TAB(iRow-1, 4);
                    end

                    %if heading, start, and time are defined, and end
                    %locations are nan, update end location.
                    if nnz(~isnan(TAB(iRow, [1,2,5,6])))==4 && nnz(isnan(TAB(iRow, [3,4])))== 2
                        TAB(iRow, 3) = ceil(...
                            TAB(iRow, 1) + 1 + ... % add extra buffer
                            cosd(TAB(iRow, 6)) * TAB(iRow, 5) * ...
                            app.AircraftVelocity.Value/app.SimulationResolutionEditField.Value...
                        );
                        app.SimulationResolutionEditField.Value
                        
                        TAB(iRow, 4) = ceil(...
                            TAB(iRow, 2) + 1+... % add extra buffer
                            sind(TAB(iRow, 6)) * TAB(iRow, 5) * ...
                            app.AircraftVelocity.Value/app.SimulationResolutionEditField.Value...
                        );
                    elseif nnz(isnan(TAB(iRow, 5:6))) == 2
                    
                        DX = TAB(iRow, 3) - TAB(iRow,1);
                        DY = TAB(iRow, 4) - TAB(iRow,2);
                        DT = ... % leg time calculation
                            sqrt(DX^2 + DY^2) * app.SimulationResolutionEditField.Value / ...
                            app.AircraftVelocity.Value ...
                        ;
                        TAB(iRow, 5) = DT;

                        HEADING = atan2d(DY, DX);
                        TAB(iRow, 6) = HEADING;
                        
                    end
                    
                    % heading calculation (math convention)
                    

                    % write to the table
                    app.UITable.Data = TAB;
                catch ME
                    ME.message
                    
                end
            end
        end
        
        function app = updateMessages(app,message, varargin)
            message = ['{',message,'}'];
            if isempty(varargin)
                app.messages{end+1} = message;
            elseif strcmp(varargin{1}, 'append')
                app.messages{end} = [app.messages{end}, message];
            end

            messagesTemp = fliplr(app.messages);
            try
                app.MessagesUI1.Text = messagesTemp(1:min(5,length(app.messages)));
                app.MessagesUI2.Text = messagesTemp(1:min(5,length(app.messages)));
                app.MessagesUI3.Text = messagesTemp(1:min(5,length(app.messages)));
            end
            
        end
        
        function [row, col] = latlonTorowcol(app, lat, lon)
            % convert from latitude/longitude to simulation idx.
          

            [~,~,d,~] = latlonTodisaz(lat, lon, app.xlat, app.xlon);
            
            [row, col] = find(d == min(d(:)));
        end

        function [lat, lon] = rowcolTolatlon(app,row, col)
            %convert from row/col in simulation space to lat/lon
           
            lat = app.xlat(round(row), round(col));
            lon = app.xlon(round(row), round(col));
        end
        
        function varLoaded = loadVar(app, varName)
            varLoaded = ncread(...
                [...            % load variable at appropriate time
                    app.FileList(app.TimeIndexEditField.Value).folder, '/', ...
                    app.FileList(app.TimeIndexEditField.Value).name ...
                ], ...
                varName...   %variable to load
            );
            
        end
        
        function app = setAltitudeSlidersAndValues(app, idx)
            if idx >= min(app.LevelToPlotSlider.Limits) && ...
                    idx <= max(app.LevelToPlotSlider.Limits)
                
            
            
                app.AltToPlotText.Value = double(app.verticalLevels.(app.AltToPlotUnits.Value)(idx));
                app.LevelToPlot.Value = idx;
                app.LevelToPlotSlider.Value = idx;

                app.UpdateFigure()
            else
                app.updateMessages(['\color{red}Selected level is outside of available levels. The range is = [', ...
                    num2str(min(app.LevelToPlotSlider.Limits)), ', ', num2str(max(app.LevelToPlotSlider.Limits)), ']'])
                app.AltToPlotText.Value = app.verticalLevels.(app.AltToPlotUnits.Value)(1);
                app.LevelToPlot.Value = 1;
                app.LevelToPlotSlider.Value = 1;
            end
        end
        
        function app = setTimeSlidersAndValues(app, idx)
            
            if idx >= min(app.TimeIndexSlider.Limits) && ...
                    idx <= max(app.TimeIndexSlider.Limits)
                app.TimeIndexSlider.Value = idx;
                app.TimeIndexEditField.Value = idx;

                % if the start time same as plotted time checkbox is
                % marked, link the two.
                if app.SameasplottedtimeCheckBox.Value == 1
                    
                    ticks = app.StartTimeSlider.MajorTicks;
                    app.StartTimeSlider.Value = ticks(idx);
                    app.StartTimeEditField.Value = ticks(idx);
                end


                % update figure
                UpdateFigure(app);
            else
                app.TimeIndexEditField.Value = 1;
            end
        end
        
        function [pathXFinal, pathYFinal] = getWarpedFlight(app, warpBool)
            if warpBool
                sTAB = size(app.UITable.Data);
                TABT = app.UITable.Data';
                pathX = TABT([1,3:sTAB(2):3+sTAB(2)*(sTAB(1)-1)]);
                pathY = TABT([2,4:sTAB(2):4+sTAB(2)*(sTAB(1)-1)]);
                pathTimes = cumsum([0, TABT(5,:)]);
    
               
                
                fileTimesShifted = app.FileTimes - app.FileTimes(app.TimeIndexSlider.Value);
    
                indEnd = find(fileTimesShifted > max(pathTimes), 1);
                if isempty(indEnd); indEnd = numel(fileTimesShifted); end
                queryTimes = fileTimesShifted(find(fileTimesShifted>=0, 1):indEnd);
                
                pathXinterp = interp1(pathTimes, pathX, queryTimes, "linear", "extrap");
                pathYinterp = interp1(pathTimes, pathY, queryTimes, "linear", "extrap");
                
                mask = ~isnan(pathYinterp);
                pathXinterp = pathXinterp(mask);
                pathYinterp = pathYinterp(mask);
                queryTimes = queryTimes(mask);
                
    
                for i=2:numel(pathXinterp)
                    if find(fileTimesShifted == queryTimes(i) ) == numel(fileTimesShifted)+1
                        pathXInterp(i) = NaN;
                        pathYinterp(i) = NaN;
                    else
                        tf = app.deformationTransforms{find(fileTimesShifted == queryTimes(i),1)-1};
                        [y,x] = transformPointsForward(tf,  pathYinterp(i:end), pathXinterp(i:end));
                        pathXinterp(i:end) = x;
                        pathYinterp(i:end) = y;
                    end
                end
                
                
                pathTimes(find(pathTimes>max(queryTimes),1))=max(queryTimes);
                pathTimes(pathTimes>max(queryTimes)) = [];
                pathXFinal = interp1(queryTimes, pathXinterp, pathTimes);
                pathYFinal = interp1(queryTimes, pathYinterp, pathTimes);
            else
                sTAB = size(app.UITable.Data);
                TABT = app.UITable.Data';
                pathXFinal = TABT([1,3:sTAB(2):3+sTAB(2)*(sTAB(1)-1)]);
                pathYFinal = TABT([2,4:sTAB(2):4+sTAB(2)*(sTAB(1)-1)]);
            end


            
        end
        
        function  writeNamelists(app)
            checkedNodes = app.Tree.CheckedNodes;
            % clear all prior namelists and scan lists
            delete("./output/namelists/*.txt")
            delete("./output/scanlists/*.txt")

            
            %% loop through all nodes, finding the parent nodes.
            for iNode = 1:numel(checkedNodes)
                if numel(checkedNodes(iNode).Parent)>0
                    try
                        %if ~exist("./output/namelists/"+lower(checkedNodes(iNode).Parent.Text(1:3)+"_namelist.txt"), 'file')
                        WriteNamelistHelper(app, checkedNodes(iNode))
                        %end
                        
                        % add relevant scan table to final scanlist.
                        fileID1 = fopen(sprintf("./.meta/scanTables/%s_%s.txt",lower(checkedNodes(iNode).Parent.Text(1:3)), lower(checkedNodes(iNode).Text)), 'r');
                        fname = sprintf("./output/scanlists/%s_%s.txt",lower(checkedNodes(iNode).Parent.Text(1:3)), lower(checkedNodes(iNode).Text));
                    
                        
                         
                        if exist(fname,'file')
                            %skip over preamble if this is the second or
                            %greater time visiting this file.
                            
                            fgetl(fileID1)
                            fgetl(fileID1)
                            
                            fgetl(fileID1)
                            fgetl(fileID1)
                            fileID2 = fopen(fname,'a');
                        else
                            fileID2 = fopen(fname,'w');
                        end
    
                        
                        
                        while true
                            line = fgetl(fileID1);
    
                            fprintf(fileID2, [line,'\n']);
                            
                            if contains(line, "</sweep>")
                                break
                            end
                        end
                        fprintf("\n")
                        fclose(fileID2);
                        fclose(fileID1);
                    catch ME
                        if isa(ME.identifier, 'MATLAB:FileIO:InvalidFid')
                            continue
                        end
                    end
                       
                end
            end

            function WriteNamelistHelper(app, node)
                % subfunction to write the text of the namelist(s)
                     
                namelistText = [...
                    "&options", ...
                    " seed = 875270836, 1299229797" ...
                ];
            
                %% get wrf_glob pattern
                try
                    globPattern = app.FileList(1).name;
                    globPattern(char(app.TimeGlobPattern.Value)~='?') = "?";
                    namelistText(end+1) = sprintf(" wrf_glob_pattern = '%s/%s'", app.FileList(1).folder,globPattern);
                catch
                    error("file search pattern is incorrect.")
                end

                %% output filename format
                try
                    formatString = app.output_filename_format_string.Value;
                    formatString = strrep(formatString, "PNL", lower(node.Parent.Text(1:3)));
                    formatString = strrep(formatString, "SCN", lower(node.Text(1:3)));
                    namelistText(end+1) = sprintf(" output_filename_format_string = %s", formatString);
                catch
                    error("output file name format is incorrect")
                end

                %% leg details
                
                TAB = app.UITable.Data;
                TABT = round(TAB');
                sTAB = size(TAB);
                try
                    % get warped flight path and add padding to the end of
                    % the path
                    [x,y] = app.getWarpedFlight(app.WarpFlightPathCheckBox.Value);
                    dist = sqrt(diff(x).^2+diff(y).^2);
                    
                    modifiedSpeed = floor(sum(dist)*app.SimulationResolutionEditField.Value/round(sum(TAB(:,5))));
                    theta = atan2d(y(end)-y(end-1), x(end)-x(end-1));
                    % add 5% of the distance of the flight path to the end
                    % of the run as padding.
                    fudgeDistance = round(sum(TAB(:,5)))*app.AircraftVelocity.Value*0.05;
                    x(end) = min(...
                        x(end) + (fudgeDistance/app.SimulationResolutionEditField.Value)*cosd(theta), ...
                        size(app.environment,1) ...
                    );
                    y(end) = min(...
                        y(end) + (fudgeDistance/app.SimulationResolutionEditField.Value)*sind(theta), ...
                        size(app.environment,2) ...
                    );
                    

                    namelistText(end + 1) = " leg_initial_time = " + app.StartTimeEditField.Value;
                    namelistText(end + 1) = " leg_time_seconds = " + string(int32(round(sum(TAB(:,5))))) + ", ";
                    
                    namelistText(end + 1) = " time_evolution = .TRUE.";
                    
                    namelistText(end + 1) = " flight_waypoints_x = " + ...
                        strjoin(...
                            string(...
                                round(x) ...
                            ), ...
                        ", "...
                    );
                    namelistText(end + 1) = " flight_waypoints_y = " + ...
                        strjoin(...
                            string(...
                                round(y) ...
                            ), ...
                        ", " ...
                    );
                    namelistText(end + 1) = " flight_waypoints_vert = " + ...
                        strjoin(...
                            string(...
                                ones(size(TABT([2,4:sTAB(2):4+sTAB(2)*(sTAB(1)-1)])))*str2double(app.flight_level_waypoints_vert.Value)...
                            ), ...
                        ", " ...
                    );

                    namelistText(end + 1) =  " flight_level_coordinate = "+app.flight_level_coordinate.Value;
                    namelistText(end) = sprintf(" air_speed = %s. ! meters per second", num2str(app.AircraftVelocity.Value));
                catch
                    app.updateMessages("\color{red}Make sure Start Time, table, and flight level coordinates are set properly")
                end

                
                namelistText(end + 1:end+4) = [...
                    " herky_jerky = .FALSE.", ...
                    " bwtype = 0,", ...
                    " ref_angle = 270,", ...
                    "/"...
                ];

                %% attitude
                namelistText(end + 1:end+5) = [...
                    "&attitude_external_source", ...
                    " use_external_attitudes = .False.,", ...
                    " attitude_file = 'attitude.nc',", ...
                    " attitude_orientation_rotate_degrees = 182.0,", ...
                    "/" ...
                ];

                %% scanning            

                namelistText(end + 1) = "&scanning";
                namelistText(end + 1) = " CRSIM_Config                  = 'CONFIG_crsim'";
                namelistText(end + 1) = sprintf(" scanning_table                = './output/scanlists/%s_%s.txt'",lower(node.Parent.Text(1:3)), lower(node.Text));
                namelistText(end + 1) = " pulse_repetition_frequency    = "+app.pulse_repetition_frequency.Value;
                namelistText(end + 1) = " pulses_per_pulse_set          = "+app.pulses_per_pulse_set.Value;
                namelistText(end + 1) = " revisits_per_acquisition_time = "+app.revisits_per_acquisition_time.Value;
                namelistText(end + 1) = " beams_per_acquisition_time    = "+app.beams_per_acquisition_time.Value;
                namelistText(end + 1) = " skip_seconds_between_scans    = "+app.skip_seconds_between_scans.Value;
                namelistText(end + 1) = " meters_between_gates          = "+app.meters_between_gates.Value;
                namelistText(end + 1) = " meters_to_center_of_first_gate= "+app.meters_to_center_of_first_gate.Value;
                namelistText(end + 1) = " max_range_in_meters           = "+app.max_range_in_meters.Value;
                namelistText(end + 1) = "/";

                %% wrap up
                namelistText(end + 1) = "&config_output";
                namelistText(end + 1) = "/";

                %% write to file
                writelines(namelistText, "./output/namelists/"+lower(node.Parent.Text(1:3)+"_"+lower(node.Text)+"_namelist.txt"))
                

            end
            
        end
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFunction(app)
            %% add helpers folder to path
            addpath(genpath("./helpers/"))
            app.updateMessages('\color{green} starting app!')

            
            %% make .meta/* or output folder if none exists
            if ~exist("./.meta", "file")
                app.updateMessages('Writing metadata folders...')
                mkdir("./.meta")
            end

            if ~exist("./.meta/scanTables", "file")
                app.updateMessages('Writing `./.meta/scanTables` folder')
                mkdir("./.meta/scanTables")
            end

            if ~exist("./output", "file")
                app.updateMessages('Writing `./output` folder')
                mkdir("./output")
            end

            if ~exist("./output/namelists", "file")
                
                mkdir("./output/namelists")
            end

            if ~exist("./output/scanlists", "file")
           
                mkdir("./output/scanlists")
            end

            

            %% save default scanlists
            app.RecalculateScanlistsValueChanged(app);
            
            
            %% try to load last state of application
            try 
                savedState = load('./.meta/savedState.mat', 'savedState');
                savedState = savedState.savedState;
                suffix = load('./.meta/savedState.mat', 'suffix');
                suffix = suffix.suffix;
                varNames = string(fieldnames(savedState));
                for iVar = 1:numel(varNames)
                    
                    
                    if ~isa(savedState.(varNames(iVar)), 'struct')
                        app.(varNames(iVar)) = savedState.(varNames(iVar));
                    else
                        subVarNames = string(fieldnames(savedState.(varNames(iVar))));
                        for jSuffix = 1:numel(subVarNames)
                            app.(varNames(iVar)).(subVarNames(jSuffix)) = savedState.(varNames{iVar}).(subVarNames(jSuffix));
                        end
                    end
                    
                end

                app.SimulationGlobPatternValueChanged(app)
                app.TimeGlobPatternValueChanged(app)
            catch ME
                

            end
            
        end

        % Value changed function: StartTimeEditField
        function StartTimeEditFieldValueChanged(app, event)
            value = app.StartTimeEditField.Value;
            [~, idx] = min(abs(value(1) - app.FileTimes));
            % if linked to plotted time, invoke that function
            if app.SameasplottedtimeCheckBox.Value == 1
                app.setTimeSlidersAndValues(idx)
            
            % Otherwise, just edit the edit field.
            else
                app.StartTimeSlider.Value = app.FileTimes(idx);
                app.StartTimeEditField.Value = app.FileTimes(idx);
            end
        end

        % Value changed function: SimulationGlobPattern
        function SimulationGlobPatternValueChanged(app, event)
            value = app.SimulationGlobPattern.Value;
          
            
            try
                
                app.FileList = dir(string(value));
                app.FileList = app.FileList(~ismember({app.FileList.name},{'.','..'}));

                if nnz(diff(strfind(app.FileList(1).name, digitsPattern)) > 1) > 1
                    selection = uiconfirm(app.UIFigure, ...
                        'It looks like the folder you selected contains a simulation files with a naming scheme not supported by AOSPRE, would you like me to try to rename them to a supported format?', ...
                        'Unsupported filenames found', ...
                        'Options', ["Yes", "No"]);
                    if strcmp(selection, "Yes")
                        try 
                            wrfoutToAOSPRE(app.FileList(1).folder, app.FileList(1).folder)
                        catch
                            % tell user to run wrfoutToAOSPRE.py manually
                            uiconfirm(app.UIFigure, ...
                                'Failed to move files automatically, manually invoke `./helpers/wrfoutToAOSPRE.py` to rename files (command has been copied to system clipboard).', 'Failed to rename files','Options',"OK")
                            clipboard("copy", ...
                                "python ./helpers/wrfoutToAOSPRE.py -in "+string(app.FileList(1).folder)+" -out "+string(app.FileList(1).folder))
                        end
                    end
                end
                
                
    
                if numel(app.FileList)>0
                    app.updateMessages(char(num2str(numel(app.FileList))+" files found."));
    
                
                    %% update variable dropdown based on file search
                    % get variables in file
                    
                    simInfo=ncinfo([app.FileList(1).folder,'/',app.FileList(1).name],"/");
                    varNames = string({simInfo.Variables.Name});
                    varDimNumber = cellfun(@nnz,{simInfo.Variables.Size});

                    app.VariableToPlot.Items = varNames(varDimNumber == 4);
                    
                    %% Set children visibility
                    SetFigureControlsVisibility(app)

                    

                    %% Try variable    
                    try
                        
                        varsToCheck = ["W", "Q", "QVAPOR"];
                        for iVar = 1:numel(varsToCheck)
                            
                            mask = strcmp(varNames, varsToCheck(iVar));
                            if nnz(mask) > 0
                                app.VariableToPlot.Value = varNames(find(mask, 1));
                                app.UpdateSlidersAndRanges()
                                app.VariableToPlotValueChanged()
                            end
                        end
                       
                    catch ME
                        ME.message
                    end

                    
                    
                    

                    %% Guess simulation spatial resolution
                    try
                        
                        XLAT = ncread(...
                            [...
                               app.FileList(1).folder, '/', ...
                               app.FileList(1).name ...
                            ], ...
                            'XLAT' ...
                        );
                        commonResolutions = [50, 100, 150, 200, 250, 300, 400, 500, 600, 700, 800, 900, 1000, 1500, 2000, 3000, 9000, 27000]; % m
                        dlat = XLAT(2,2,1)-XLAT(1,1,1);
                        dy = dlat*111131.745;
                       
                        [~,bestGuessResolutionIdx] = min(abs(dy-commonResolutions));
                        bestGuessResolution = commonResolutions(bestGuessResolutionIdx);
                        app.SimulationResolutionEditField.Value = bestGuessResolution;
                        app.updateMessages(['Simulation resolution guess: ', num2str(bestGuessResolution), ' m'])
                    end

                    %% Guess time filename pattern
                    try
                        app.TimeGlobPatternValueChanged(app)
                    catch ME
                        ME.message
                    end

                    %% Update figure MMMM d, yyyy HH:mm:ss
                    UpdateFigure(app)
                    
                    %% move to next panel if everything passes
                    if app.SimulationResolutionEditField.Value ~= 1
                        app.updateMessages('\color{green} All initial checks passed, moving to flight planning tab!')
                        app.TabGroup.SelectedTab = app.FlightPlanningTab;
                    end
                else
                    app.updateMessages('\color {orange} Choose directory with at least 2 files.')
                end
            catch ME
                ME.message
                app.updateMessages({'\color{red} Search Failed, try again.'});
            end
            


            function wrfoutToAOSPRE(inputDir, outputDir)
            
                % Get files in wrfoutDir
                wrfFiles = dir(fullfile(inputDir, 'wrfout*'));
                wrfFiles = {wrfFiles.name};
                
                % Sort wrfFiles
                wrfFiles = sort(wrfFiles);
            
                % Find start time
                startString = split(wrfFiles{1}, '_');
                timeString = strcat(startString{3}, '_', startString{4});
                startDt = datetime(timeString, 'InputFormat', 'yyyy-MM-dd_HH:mm:ss');
                
                % Convert to seconds since start time for each file and rename
                for i = 1:length(wrfFiles)
                    file = wrfFiles{i};
                    if ~contains(file, 'wrfout')
                        continue;
                    end
                    
                    fileString = split(file, '_');
                    prefix = strcat(fileString{1}, '_', fileString{2}, '_');
                    timeString = strcat(fileString{3}, '_', fileString{4});
                    
                    % Check if file is wrfout file
                    fileTime = datetime(timeString, 'InputFormat', 'yyyy-MM-dd_HH:mm:ss');
                    
                    % Calculate time difference
                    timeDelta = fileTime - startDt;
                    secondsdiff = seconds(timeDelta);
                    
                    % Rename file
                    newFile = sprintf('%s%06.0f.nc', prefix, secondsdiff);
                    sourceFile = fullfile(inputDir, file);
                    destFile = fullfile(outputDir, newFile);
                    
                    command = sprintf('copyfile(''%s'', ''%s'')', sourceFile, destFile);
                    fprintf(command);
                    movefile(sourceFile, destFile);
                end
            end
        end

        % Value changed function: VariableToPlot
        function VariableToPlotValueChanged(app, event)
            value = app.VariableToPlot.Value;
            % update sliders
            UpdateSlidersAndRanges(app)
            % update figure
            UpdateFigure(app);
       
        end

        % Value changed function: TimeIndexSlider
        function TimeIndexSliderValueChanged(app, event)
            % link edit field to slider value
            value = app.TimeIndexSlider.Value;

            app.setTimeSlidersAndValues(round(value));
            

            % update figure
            UpdateFigure(app);
        end

        % Value changed function: TimeIndexEditField
        function TimeIndexEditFieldValueChanged(app, event)
            % link slider value to edit field
            value = app.TimeIndexEditField.Value;
            
            app.setTimeSlidersAndValues(value)
            
        end

        % Value changed function: LevelToPlotSlider
        function LevelToPlotSliderValueChanged(app, event)
            %link edit field to slider value
            value = app.LevelToPlotSlider.Value;
            app.LevelToPlotSlider.Value = round(value);
            
            app.setAltitudeSlidersAndValues(app.LevelToPlotSlider.Value);

            % update figure
            UpdateFigure(app);
        end

        % Value changed function: LevelToPlot
        function LevelToPlotValueChanged(app, event)
            value = app.LevelToPlot.Value;
            app.LevelToPlot.Value = round(value);

            app.setAltitudeSlidersAndValues(app.LevelToPlot.Value)

                
        end

        % Button pushed function: ManuallyaddflightlegButton
        function ExtractWaypoints(app, event)
            % update table if new marks are added to map

            if isempty(app.UITable.Data) | (nnz(isnan(app.UITable.Data)) == numel(app.UITable.Data))
                app.flightPath = drawpolyline(app.UIAxes);
                app.flightPath.addlistener('ROIMoved', @(~,~) app.UpdateTableLogic())
            else
                x = app.UITable.Data(:,[1,3]);
                x = x(:);
                y = app.UITable.Data(:,[2,4]);
                y = y(:);

                xc = x(~isnan(x) | ~isnan(y));
                yc = y(~isnan(x) | ~isnan(y)); % simulation space

                [yf, xf] = app.rowcolTolatlon(xc(end), yc(end)); %simulation space to lat/lon
                
                app.flightPath = images.roi.Polyline(app.UIAxes);
                beginDrawingFromPoint(app.flightPath, [xf, yf])
            end
            
            %convert flightPath back to simulation space
            for ic = 1:size(app.flightPath.Position,1)
                lon = app.flightPath.Position(ic,1);
                lat = app.flightPath.Position(ic,2);
                [x,y] = app.latlonTorowcol(lat, lon);
                app.flightPath.Position(ic,:) = [x,y];
            end
            
            app.UpdateTableLogic()
            app.UpdateFigure()
            
            
    
       
        end

        % Value changed function: TimeGlobPattern
        function TimeGlobPatternValueChanged(app, event)
            value = app.TimeGlobPattern.Value;
            %% try to guess filenaming scheme (if defaults aren't written)
            if ~contains(value,"?")
                datetimeChars = 'yyMMddHHmmss';
                
                COUNTER = 1;
                timePatternGuess = app.FileList(1).name;
                %find the columns where fileNames change
                
                
                inds = strfind(timePatternGuess, digitsPattern);
                inds = inds(end-5:end);
                allInds = 1:length(timePatternGuess);

                app.TimeGlobPattern.Value = timePatternGuess;
                app.TimeGlobPattern.Value(~ismember(allInds,inds)) = '?';
                app.TimeGlobPattern.Value(allInds(inds)) = 's';
                
               
            end

            %% update filetimes variable
            datetimeChars = 'yMdHms';
            
            findNumber = @(jChar, iFile) str2num(app.FileList(iFile).name(app.TimeGlobPattern.Value == datetimeChars(jChar)));
            
            for iFile = 1:numel(app.FileList)
                DT = 0;
                for jChar = 1:numel(datetimeChars)
                    
                    dur = findNumber(jChar,iFile);
                    if ~isempty(dur)
                        DT = DT + findNumber(jChar,iFile);
                    end
                end
                app.FileTimes(iFile) = DT;
            end

            %% Update the start time slider
            app.StartTimeSlider.Limits = [min(app.FileTimes), max(app.FileTimes)];
            app.StartTimeSlider.MajorTicks = app.FileTimes;
            app.StartTimeSlider.MinorTicks = [];
            app.StartTimeSlider.Value = min(app.FileTimes);
            app.StartTimeEditField.Value = min(app.FileTimes);
            
        end

        % Callback function
        function UITableCellEdit(app, event)
            indices = event.Indices;
            newData = event.NewData;

            %% try to complete flight calculation for each row

            app.UpdateTableLogic()
            app.UpdateFigure()
        end

        % Value changed function: SimulationResolutionEditField
        function SimulationResolutionEditFieldValueChanged(app, event)
            value = app.SimulationResolutionEditField.Value;
            app.UpdateTableLogic()

            if ~isempty(app.SimulationGlobPattern.Value)
                app.TabGroup.SelectedTab = app.FlightPlanningTab;
            end
        end

        % Value changed function: StartTimeSlider
        function StartTimeSliderValueChanged(app, event)
            value = app.StartTimeSlider.Value;
            [~, idx] = min(abs(value(1) - app.FileTimes));
            % if linked to plotted time, invoke that function
            if app.SameasplottedtimeCheckBox.Value == 1
                app.setTimeSlidersAndValues(idx)
            
            % Otherwise, just edit the edit field.
            else
                app.StartTimeSlider.Value = app.FileTimes(idx);
                app.StartTimeEditField.Value = app.FileTimes(idx);
            end
            
            
        end

        % Menu selected function: ClearSelection
        function ClearSelectionSelected(app, event)
            %clear selection
            
            for iCell = 1:size(app.UITable.Selection,1)
                app.UITable.Data(...
                    app.UITable.Selection(iCell, 1), ...
                    app.UITable.Selection(iCell, 2) ...
                ) = NaN;
            end
            
            %% update logic once cleared
            app.UpdateTableLogic()
            app.UpdateFigure()

        end

        % Menu selected function: AddnewLegMenu
        function AddnewLegMenuSelected(app, event)
            % add a new leg to the table
            app.UITable.Data(end+1, :) = NaN([1,6]);
        end

        % Menu selected function: LegRowMenu
        function LegRowMenuSelected(app, event)
            
            app.UITable.Data(...
                    app.UITable.Selection(:,1),:) = [] ...
            ;

            app.UpdateFigure()

        end

        % Menu selected function: RefreshSelectionMenu
        function RefreshSelectionMenuSelected(app, event)
            %% refresh entire row if start is selected
            % delete pairs of cells
            
            pairs = {[],[],  [3,4], [3,4], [5,6], [5,6]};
            for iCell = 1:size(app.UITable.Selection,1)
                app.UITable.Data(...
                    app.UITable.Selection(iCell, 1), ...
                    pairs{app.UITable.Selection(iCell, 2)} ...
                ) = NaN;
            end

            %% update table after removing data
            app.UpdateTableLogic()
            app.UpdateFigure()
            
        end

        % Value changed function: AftEditField, 
        % ...and 3 other components
        function RecalculateScanlistsValueChanged(app, event)
            % family of functions to generate scanlists based on
            % defaults/inputs and save them to text files.
            app.updateMessages('Recalculating scanlists... ')
            value = app.AngularResolutionEditField.Value;

            % call the functions
            GenerateScanlists(app)
            SaveScanlists(app)
            app.updateMessages('done.', 'append')

            function [] = GenerateScanlists(app)
                SC = app.ScanTables;    % temp variable for scan tables
                                
                %% Set Nadir direction of each panel (Cartesian; y=1 is front)
                SC.rhs.origin = [1,0,0]; 
                SC.lhs.origin = [-1,0,0]; 
                SC.top.origin = [0,0,1];
                SC.bot.origin = [0,-cosd(70),-sind(70)];

                app.PanelLocations = SC;

                %% lhs/rhs 
                % rhi
                rhiTilt = -app.MaximumTiltRotEditField.Value:app.AngularResolutionEditField.Value:app.MaximumTiltRotEditField.Value;
                panelNames = ["rhs", "lhs"];
                for iPanel = 1:2
                    SC.(panelNames(iPanel)).rhi.tilt = ...
                        atan2d(...
                            SC.(panelNames(iPanel)).origin(3), ...
                            sqrt(SC.(panelNames(iPanel)).origin(1)^2 + SC.(panelNames(iPanel)).origin(2)^2)...
                        ) + rhiTilt;
                    SC.(panelNames(iPanel)).rhi.rot = ...
                        atan2d(...
                            SC.(panelNames(iPanel)).origin(2), ...
                            SC.(panelNames(iPanel)).origin(1)...
                        )*ones(size(SC.(panelNames(iPanel)).rhi.tilt));
                    
                end
    
                % ppi
                ppiTilt = -app.MaximumTiltRotEditField.Value:app.AngularResolutionEditField.Value:app.MaximumTiltRotEditField.Value;
                for iPanel = 1:2
                    SC.(panelNames(iPanel)).ppi.rot = ...
                        atan2d(...
                            SC.(panelNames(iPanel)).origin(2), ...
                            SC.(panelNames(iPanel)).origin(1)...
                        ) + ppiTilt;
                    SC.(panelNames(iPanel)).ppi.tilt = ...
                        atan2d(...
                            SC.(panelNames(iPanel)).origin(3), ...
                            sqrt(SC.(panelNames(iPanel)).origin(1)^2 + SC.(panelNames(iPanel)).origin(2)^2)...
                        )*ones(size(SC.(panelNames(iPanel)).ppi.rot));
                end
                % fore
                
                for iPanel = 1:2
                    SC.(panelNames(iPanel)).fore.rot = ...
                        SC.(panelNames(iPanel)).rhi.rot + ...
                        app.ForeEditField.Value *sign(SC.(panelNames(iPanel)).origin(1)); %flip sign if on lhs vs rhs.
                    SC.(panelNames(iPanel)).fore.tilt = SC.(panelNames(iPanel)).rhi.tilt;
                end
               
    
                % aft
                for iPanel = 1:2
                    SC.(panelNames(iPanel)).fore.rot = ...
                        SC.(panelNames(iPanel)).rhi.rot - ...
                        app.ForeEditField.Value * sign(SC.(panelNames(iPanel)).origin(1)); %flip sign if on lhs vs rhs.
                    SC.(panelNames(iPanel)).fore.tilt = SC.(panelNames(iPanel)).rhi.tilt;
                end
                
                %% top/bot
                % fore
                ppiTilt = 90+(-app.MaximumTiltRotEditField.Value:app.AngularResolutionEditField.Value:app.MaximumTiltRotEditField.Value);
                panelNames = ["top", "bot"];
                for iPanel = 1:2
                    
                    SC.(panelNames(iPanel)).fore.rot = ...
                        atan2d(...
                            SC.(panelNames(iPanel)).origin(2), ...
                            SC.(panelNames(iPanel)).origin(1)...
                        ) + ppiTilt;
    
                    SC.(panelNames(iPanel)).fore.tilt = ...
                        atan2d(...
                            SC.(panelNames(iPanel)).origin(3), ...
                            SC.(panelNames(iPanel)).origin(2) ...
                        ) - ... 
                        app.ForeEditField.Value * sign(SC.(panelNames(iPanel)).origin(3))*ones(size(SC.(panelNames(iPanel)).fore.rot)); %flip sign if on top vs. bot
                end
                
                % aft
                panelNames = ["top", "bot"];
                ppiTilt = 270+(-app.MaximumTiltRotEditField.Value:app.AngularResolutionEditField.Value:app.MaximumTiltRotEditField.Value);
                for iPanel = 1:2
                    SC.(panelNames(iPanel)).aft.rot = ...
                        atan2d(...
                            SC.(panelNames(iPanel)).origin(2), ...
                            SC.(panelNames(iPanel)).origin(1)...
                        ) + ppiTilt;
    
                    SC.(panelNames(iPanel)).aft.tilt = ...
                        atan2d(...
                            SC.(panelNames(iPanel)).origin(3), ...
                            SC.(panelNames(iPanel)).origin(2) ...
                        ) - ... 
                        app.AftEditField.Value * sign(SC.(panelNames(iPanel)).origin(3))*ones(size(SC.(panelNames(iPanel)).aft.rot)); %flip sign if on top vs. bot
                end
    
    
                %% convert to aircraft coordinates
                % rotation for lhs and rhs are meteo convention, tilt is fine.
                panelNames = ["lhs", "rhs", "bot","top"];
                for iPanel = 1:numel(panelNames)
                    scanNames = string(fieldnames(SC.(panelNames(iPanel))));
                    scanNames(contains(scanNames, 'origin')) = [];
                    
                    
                    for jScans = 1:numel(scanNames)
                        SC.(panelNames(iPanel)).(scanNames(jScans)).rot = wrapTo360(90-SC.(panelNames(iPanel)).(scanNames(jScans)).rot);
                    end
                end
    
                %% save temp variable back into app property
                app.ScanTables = SC;

            
            end
            function [] = SaveScanlists(app)
                % generate scan tables from app.ScanTables and save into text
                % file
                %% Store scanlists in text files using Fortran I/O format.
                SC = app.ScanTables;
                preamble = [...
                    "# see https://github.com/NCAR/AOSPRE/tree/main/Scanning_Table_Library for more information", ...
                    "PRIMARY_AXIS = Z" ...
                ];
    
                panelNames = ["lhs", "rhs","top","bot"];
                for iPanel = 1:numel(panelNames)
                    scanNames = string(fieldnames(SC.(panelNames(iPanel))));
                    scanNames(contains(scanNames, 'origin')) = [];
                    for jScan = 1:numel(scanNames)
                        %% generate filename
                        
                        scanlistFileName = "./.meta/scanTables/"; %init filename
                        scanlistFileName = scanlistFileName +panelNames(iPanel)+"_"+scanNames(jScan)+".txt";
    
                        %% generate text to print
                        textToPrint = preamble;
    
                        textToPrint(end+1) = "SWEEP_MODE = "+upper(scanNames(jScan));
                        textToPrint(end) = strrep(textToPrint(end), 'PPI', 'sector');
                        textToPrint(end+1) = "PARAMETERS = ROT , TILT";
                
                        
                        textToPrint(end + 1) = "<sweep>";
                        for kBeams = 1:numel(SC.(panelNames(iPanel)).(scanNames(jScan)).rot)
                            rot = SC.(panelNames(iPanel)).(scanNames(jScan)).rot(kBeams);
                            tilt = SC.(panelNames(iPanel)).(scanNames(jScan)).tilt(kBeams);
                            textToPrint(end+1) = string(sprintf(" %3.1f %3.1f",rot,tilt));
                        end
                        textToPrint(end+1) = "</sweep>";
                        writelines(textToPrint, scanlistFileName)
    
                    end
    
                end
    
    
                
            end
            
        end

        % Callback function: Tree
        function GenerateNamelists(app, event)
            app.writeNamelists()
            
            
        end

        % Button pushed function: InitAOSPREButton
        function InitAOSPREButtonPushed(app, event)
            
            %% Check to make sure everything has been filled in. 
            
            varsToCheck = {...
                'TimeGlobPattern', ...
                'AircraftVelocity', ...
                'TimeGlobPattern', ...
                'SimulationGlobPattern', ...
                'AftEditField', ...
                'ForeEditField', ...
                'MaximumTiltRotEditField', ...
                'AngularResolutionEditField', ...
                'StartTimeEditField', ...
                'max_range_in_meters', ...
                'meters_to_center_of_first_gate', ...
                'beams_per_acquisition_time', ...
                'beams_per_acquisition_time', ...
                'revisits_per_acquisition_time', ...
                'pulses_per_pulse_set', ...
                'pulse_repetition_frequency', ...
                'flight_level_waypoints_vert', ...
                'flight_level_coordinate', ...
                'output_filename_format_string'
            };

            pathsToCheck = {...
                'AOSPREPath', ...
                'CRSIM_ConfigPath'
            };
            errorCounter = 0;
            for iVar = 1:numel(varsToCheck)
                if isempty(app.(varsToCheck{iVar}).Value)
                    % find all parents
                    findParentsHelper = app.(varsToCheck{iVar});
                    parentsString = [];
                    while ~isa(findParentsHelper, 'matlab.ui.container.Tab')
                        findParentsHelper = findParentsHelper.Parent;
                        parentsString = [findParentsHelper.Title, '/', parentsString];
                    end
                    errorCounter = errorCounter + 1;
                    app.updateMessages(['\color{red} Variable "', parentsString, varsToCheck{iVar}, '" is missing.'])
                    
                    
                end
            end
            
            for iPath = 1:numel(pathsToCheck)
                if ~isa(app.(pathsToCheck{iPath}).Value, 'char')
                    % find all parents
                    findParentsHelper = app.(pathsToCheck{iPath});
                    parentsString = [];
                    while ~isa(findParentsHelper, 'matlab.ui.container.Tab')
                        findParentsHelper = findParentsHelper.Parent;
                        parentsString = [findParentsHelper.Title, '/', parentsString];
                    end

                    varsToCheck{iPath}.Parent
                    app.updateMessages(['\color{red} Variable "', parentsString, app.(pathsToCheck{iPath}).Title, '" is missing.']);
                    errorCounter = errorCounter + 1;
                end
            end

            

            if errorCounter > 0
                return
            end

            %% write namelists
            app.writeNamelists()   

            
            %% write bash script
            textToWrite = "ulimit -s unlimited";
            textToWrite(end+1) = sprintf("AOSPATH='%s'",app.AOSPREPath.Value);


            namelists = dir("./output/namelists/*.txt");
            namelistText = "";
            for i = 1:numel(namelists)
                namelistText = namelistText+[namelists(i).folder, '/',namelists(i).name, ' '];
            end
            textToWrite(end+1) = sprintf("NAMELISTS='%s'", namelistText);
            textToWrite(end+1) = "$AOSPATH $NAMELISTS > ./output/aospre.log 2>&1 &";
            writelines(textToWrite, "./output/run-aospre.sh")
            
            app.updateMessages('\color{green} AOSPRE-GUI completed without errors!')
            app.updateMessages('Run `bash ./output/run-aospre.sh` in terminal to run AOSPRE executable (command copied to clipboard).')
            clipboard('copy','bash ./output/run-aospre.sh')


        end

        % Close request function: UIFigure
        function shutdownFunction(app, event)


            variablesToSave = {   "SimulationGlobPattern",    "TimeGlobPattern",  "SimulationResolutionEditField",  "AOSPREPath", "UITable", "deformationTransforms", "Targets", "WarpTargetsCheckBox","WarpFlightPathCheckBox", "skip_seconds_between_scans", "max_range_in_meters", "meters_to_center_of_first_gate", "meters_between_gates", "CRSIM_ConfigPath", "beams_per_acquisition_time","revisits_per_acquisition_time", "pulses_per_pulse_set","pulse_repetition_frequency"};
            suffix = {                      "Value",                    "Value",            "Value",                          "Value",      "Data",    "",            "",        ["Value", "Visible"], ["Value", "Visible"], "Value", "Value", "Value", "Value", "Value", "Value","Value", "Value","Value"};
            
            for iVar = 1:numel(variablesToSave)
                if any(suffix{iVar} == "")
                    savedState.(variablesToSave{iVar}) = app.(variablesToSave{iVar});
                else
                    for jSuffix = 1:numel(suffix{iVar})
                        savedState.(variablesToSave{iVar}).(suffix{iVar}(jSuffix)) = app.(variablesToSave{iVar}).(suffix{iVar}(jSuffix));
                    end
                end

            end


            save('./.meta/savedState.mat', "savedState", "suffix")

            delete(app)
        end

        % Button pushed function: InteractiveTargetPlannerButton
        function InteractiveTargetPlannerButtonPushed(app, event)
            %% call setTargetInteractive class
            app.updateMessages('Run interactive target acquisition app.')
            TI = setTargetInteractive(app.environment);
            uiwait(TI.fig)
            
            %% save target(s) information
            if ~isempty(TI.Targets.map)
                
                if isempty(app.Targets.map)
                    app.Targets.map = 0;
                end
                app.Targets.map = max(app.Targets.map(:)) + TI.Targets.map;
                app.Targets.weight = [app.Targets.weight, TI.Targets.weight];
                app.Targets.dwellTime = [app.Targets.dwellTime, TI.Targets.dwellTime];
                app.Targets.loadedBool = 1;
                % ensure that the targets are only labeled with non-zero integers.
                app.Targets.map = app.Targets.map - min(app.Targets.map);


                

                
                %save('./.meta/targets.mat', 'Targets')
                
                
                %% update figure once targets have been retrieved
                app.UpdateFigure()
            end
            
        end

        % Menu selected function: RefreshMenu_2
        function RefreshMenu_2Selected(app, event)
            app.UpdateFigure()
        end

        % Button pushed function: InteractiveFlightPlannerButton
        function InteractiveFlightPlannerButtonPushed(app, event)
            
            %% initialize FlightPlanner
            FP = FlightPlanner;

            %% get environment information
            w = app.loadVar('W');
            
            wSlice = w(:,:,app.LevelToPlot.Value)';
            
            FP = initializeEnvironment(FP, wSlice, app.SimulationResolutionEditField.Value);
            if strcmp(app.StartOrContinueFlight.Value, 'Continue') | strcmp(app.StartOrContinueFlight.Value, 'Optimize')
                % pass starting location to FlightPlanner
                x0 = app.UITable.Data(:,[1,3]);
                x = [x0(1,1), x0(:,2)'];
                y0 = app.UITable.Data(:,[2,4]);
                y = [y0(1,1), y0(:,2)'];
                
                FP.priorPath = [x',y'];
            else
                bool = 1;
                if ~isempty(app.UITable.Data)
                    uic = uiconfirm(app.UIFigure, 'Are you sure? This will delete your current flight path.', 'Delete flight path');
                    bool = strcmp(uic, 'OK');
                end

                if bool == 1
                    app.UITable.Data = [];
                else
                    return
                end

                
            end
            
            FP.TrainingOptions.jobType = app.StartOrContinueFlight.Value;

            %% get aircraft information
            % aircraft speed
            FP.Aircraft.speed = app.AircraftVelocity.Value;

            % radar locations
            panelNames = ["lhs", "rhs","top","bot"];
            panelsUsed = strings(0);
            checkedNodes = app.Tree.CheckedNodes;

            FP.searchLength = max(120*120,5*app.SimulationResolutionEditField.Value);

            for iNode = 1:numel(checkedNodes)
                try
                    ind = find(panelNames == lower(checkedNodes(iNode).Parent.Text(1:3)));

                    panelsUsed(end+1) = panelNames(ind);
                catch ME
                    
                end
            end
            
            if ~isempty(panelsUsed)
                Radars = struct();
                
                for i = 1:numel(panelsUsed)
                    panelOrientation = app.ScanTables.(panelsUsed(i)).origin;
                    Radars.(panelsUsed(i)) = atan2d(-1*panelOrientation(1), panelOrientation(2));
                end
                FP.Aircraft.Radars = Radars;
            end
            
            FP.Targets = app.Targets;
            FP.Targets.map = FP.Targets.map-min(FP.Targets.map);

            FP.drawGui();

            
            uiwait(FP.fig)

            if strcmp(app.StartOrContinueFlight.Value, 'Optimize')
                app.UITable.Data = [];
            end

            app.flightPath = FP.flightPath.roi;
            app.flightPath.Parent = app.UIAxes;
            app.UpdateTableLogic()
            app.UpdateFigure();


        end

        % Button pushed function: EditflightlegverticesButton
        function EditflightlegverticesButtonPushed(app, event)
            
            if ~strcmp(app.EditflightlegverticesButton.Text, 'Done editing')
                x0 = app.UITable.Data(:,[1,3]);
                x = [x0(1,1), x0(:,2)'];
                y0 = app.UITable.Data(:,[2,4]);
                y = [y0(1,1), y0(:,2)'];
                
                
                % transform from row/col to lat/lon
                for i=1:numel(x)
                    [yt,xt] = app.rowcolTolatlon(x(i), y(i));
                    x(i) = xt;
                    y(i) = yt;
                end
                
                app.flightPath = images.roi.Polyline;
                app.flightPath.Position = [x',y'];
                app.flightPath.Parent = app.UIAxes;
                app.EditflightlegverticesButton.Text = 'Done editing';
                app.EditflightlegverticesButton.BackgroundColor = [0.5, 1, 0.5];
            else
                app.EditflightlegverticesButton.Text = 'Edit path';
                app.EditflightlegverticesButton.BackgroundColor = [0.95, 0.95, 0.95];

                %transform from lat/lon to simultion space.
                for ic = 1:size(app.flightPath.Position,1)
                    lon = app.flightPath.Position(ic,1);
                    lat = app.flightPath.Position(ic,2);
                    [x,y] = app.latlonTorowcol(lat, lon);
                    app.flightPath.Position(ic,:) = [x,y];
                end
                app.UITable.Data = [];
                app.UpdateTableLogic()
                app.UpdateFigure()
            end


            
        end

        % Button pushed function: FileBrowserButton
        function FileBrowserButtonPushed(app, event)
            dInit = char(app.SimulationGlobPattern.Value);
            if length(dInit > 1)
                
                slashInds = find(dInit == '/');
                dInit = dInit(1:slashInds(end));
            else
                dInit = './';
            end
            d = uigetdir(dInit, 'Open directory containing simulation NetCDF4 files.');
            app.SimulationGlobPattern.Value = [d,'/wrfout*'];
            
        end

        % Value changed function: AltToPlotText
        function AltToPlotTextValueChanged(app, event)
            value = app.AltToPlotText.Value;
            [~,idx] = min(abs(app.verticalLevels.(app.AltToPlotUnits.Value) - value));
            app.setAltitudeSlidersAndValues(idx);
        end

        % Button pushed function: CalculateBackgroundMotionButton
        function CalculateDeformationButtonPushed(app, event)
            selection = 'N/A';
            if exist('./.meta/deformationTransforms.mat', 'file')
                savedTransforms = load('./.meta/deformationTransforms.mat');
                selection = uiconfirm(app.UIFigure, ...
                    "Transformations have been calculated for variable '"+savedTransforms.variableCalculatedOn+"' on "+string(savedTransforms.timeOfCalculation)+". Load it instead of recalculating?","Load Saved Transforms?", ...
                    "Icon","warning","Options",["OK", "No, recalculate transformations."], "CancelOption",2);
            end

            if strcmp(selection, 'OK')
                app.deformationTransforms = savedTransforms.transforms;
                app.WarpTargetsCheckBox.Visible = 'on';
                app.WarpFlightPathCheckBox.Visible = 'on';
            else
                app.deformationTransforms = {};
                app.deformationOriginTimeStep = app.TimeIndexEditField.Value;
                % loop through each timestep and calculate the
                % registration between each.
                f = waitbar(0, 'Calculating deformations between time steps...');
                %keyboard
                for i = 1:numel(app.FileList) - 1
                    waitbar(i/(numel(app.FileList)-1), f)
                    w0 = ncread( ...
                        [app.FileList(i).folder,'/', app.FileList(i).name], ...
                        app.VariableToPlot.Value, ...
                        [1,1,app.LevelToPlot.Value,1], [inf,inf,1,1], [1,1,1,1]...
                    );
                    w1 = ncread( ...
                        [app.FileList(i+1).folder,'/', app.FileList(i+1).name], ...
                        app.VariableToPlot.Value, ...
                        [1,1,app.LevelToPlot.Value,1], [inf,inf,1,1], [1,1,1,1]...
                    );
                
                    %keyboard
                    %
                    %figure
                    %[optimizer,metric] = imregconfig("multimodal");
                    %optimizer.InitialRadius = optimizer.InitialRadius*10;
                    %tf = imregtform(normalize(w0)*255,normalize(w1)*255,"affine",optimizer,metric)
                    %
                    


                    %D=imregdemons(normalize(w0)*255, normalize(w1)*255, 15, "DisplayWaitbar",false);
                    %[x,y] = meshgrid(1:size(w0,1), 1:size(w1,2));
                    %x2 = x+D(:,:,2)';
                    %y2 = y+D(:,:,1)';

                    %mask = (w0>prctile(w0,90));
                    %tf = fitgeotform2d([x2(mask), y2(mask)], [x(mask),y(mask)], "similarity")
                    %w01=imwarp(w0, tf, 'nearest', 'OutputView',imref2d(size(w0)));

                    %imshowpair(w01,w1)


                    %optimizer.MaximumIterations = 300;
                    %tf = imregcorr(normalize(w0)*255, normalize(w1)*255, 'rigid');
                    w02 = w0;
                    w02(w02<prctile(w02,90))=0;
                    w12 = w1;
                    w12(w12<prctile(w12,90))=0;
                    tf = imregcorr(w02, w12, 'rigid',  'Window',1);

                    
                    app.deformationTransforms{end+1} = tf;
                    %app.deformationTransforms{end+1} = D;
                end
                delete(f)
    
                % set visibility of warping checkboxes only after successful
                % transform calculation.
                app.WarpTargetsCheckBox.Visible = 'on';
                app.WarpFlightPathCheckBox.Visible = 'on';
    
    
                timeOfCalculation = datetime('now');
                transforms = app.deformationTransforms;
                variableCalculatedOn = app.VariableToPlot.Value;
    
                save('./.meta/deformationTransforms', "timeOfCalculation", "transforms", "variableCalculatedOn")
            end
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 922 516];
            app.UIFigure.Name = 'MATLAB App';
            app.UIFigure.CloseRequestFcn = createCallbackFcn(app, @shutdownFunction, true);

            % Create UIAxes
            app.UIAxes = uiaxes(app.UIFigure);
            title(app.UIAxes, 'Title')
            xlabel(app.UIAxes, 'Longitude')
            ylabel(app.UIAxes, 'Latitude')
            zlabel(app.UIAxes, 'Z')
            app.UIAxes.Position = [500 200 422 317];

            % Create PlotDetailsPanel
            app.PlotDetailsPanel = uipanel(app.UIFigure);
            app.PlotDetailsPanel.Title = 'Plot Details';
            app.PlotDetailsPanel.Position = [499 0 424 198];

            % Create GridLayout3
            app.GridLayout3 = uigridlayout(app.PlotDetailsPanel);
            app.GridLayout3.ColumnWidth = {85, 40, 45, 120, 80};
            app.GridLayout3.RowHeight = {22, 22, 22, 22, 38};
            app.GridLayout3.ColumnSpacing = 7.42857142857143;
            app.GridLayout3.RowSpacing = 5.83333333333333;
            app.GridLayout3.Padding = [7.42857142857143 5.83333333333333 7.42857142857143 5.83333333333333];

            % Create TimeIndexEditField
            app.TimeIndexEditField = uieditfield(app.GridLayout3, 'numeric');
            app.TimeIndexEditField.ValueChangedFcn = createCallbackFcn(app, @TimeIndexEditFieldValueChanged, true);
            app.TimeIndexEditField.Visible = 'off';
            app.TimeIndexEditField.Layout.Row = 4;
            app.TimeIndexEditField.Layout.Column = 2;
            app.TimeIndexEditField.Value = 1;

            % Create TimeIndexLabel
            app.TimeIndexLabel = uilabel(app.GridLayout3);
            app.TimeIndexLabel.HorizontalAlignment = 'right';
            app.TimeIndexLabel.Layout.Row = 4;
            app.TimeIndexLabel.Layout.Column = 1;
            app.TimeIndexLabel.Text = 'Time Index';

            % Create TimeIndexSlider
            app.TimeIndexSlider = uislider(app.GridLayout3);
            app.TimeIndexSlider.ValueChangedFcn = createCallbackFcn(app, @TimeIndexSliderValueChanged, true);
            app.TimeIndexSlider.Visible = 'off';
            app.TimeIndexSlider.Layout.Row = 5;
            app.TimeIndexSlider.Layout.Column = [1 4];
            app.TimeIndexSlider.Value = 59.6540178571429;

            % Create LevelToPlotSlider
            app.LevelToPlotSlider = uislider(app.GridLayout3);
            app.LevelToPlotSlider.Orientation = 'vertical';
            app.LevelToPlotSlider.ValueChangedFcn = createCallbackFcn(app, @LevelToPlotSliderValueChanged, true);
            app.LevelToPlotSlider.Visible = 'off';
            app.LevelToPlotSlider.Layout.Row = [1 5];
            app.LevelToPlotSlider.Layout.Column = 5;

            % Create LevelToPlot
            app.LevelToPlot = uieditfield(app.GridLayout3, 'numeric');
            app.LevelToPlot.ValueChangedFcn = createCallbackFcn(app, @LevelToPlotValueChanged, true);
            app.LevelToPlot.Visible = 'off';
            app.LevelToPlot.Layout.Row = 3;
            app.LevelToPlot.Layout.Column = 2;
            app.LevelToPlot.Value = 1;

            % Create SimLevelLabel
            app.SimLevelLabel = uilabel(app.GridLayout3);
            app.SimLevelLabel.HorizontalAlignment = 'right';
            app.SimLevelLabel.Layout.Row = 3;
            app.SimLevelLabel.Layout.Column = 1;
            app.SimLevelLabel.Text = 'Sim. Level';

            % Create PlottedVariableDropDownLabel
            app.PlottedVariableDropDownLabel = uilabel(app.GridLayout3);
            app.PlottedVariableDropDownLabel.HorizontalAlignment = 'right';
            app.PlottedVariableDropDownLabel.Layout.Row = 1;
            app.PlottedVariableDropDownLabel.Layout.Column = 1;
            app.PlottedVariableDropDownLabel.Text = 'Plotted Variable';

            % Create VariableToPlot
            app.VariableToPlot = uidropdown(app.GridLayout3);
            app.VariableToPlot.Items = {'Set Search Pattern'};
            app.VariableToPlot.ValueChangedFcn = createCallbackFcn(app, @VariableToPlotValueChanged, true);
            app.VariableToPlot.Visible = 'off';
            app.VariableToPlot.Layout.Row = 1;
            app.VariableToPlot.Layout.Column = 2;
            app.VariableToPlot.Value = 'Set Search Pattern';

            % Create AltToPlotText
            app.AltToPlotText = uieditfield(app.GridLayout3, 'numeric');
            app.AltToPlotText.ValueChangedFcn = createCallbackFcn(app, @AltToPlotTextValueChanged, true);
            app.AltToPlotText.Visible = 'off';
            app.AltToPlotText.Layout.Row = 2;
            app.AltToPlotText.Layout.Column = 2;
            app.AltToPlotText.Value = 1;

            % Create SimLevelLabel_2
            app.SimLevelLabel_2 = uilabel(app.GridLayout3);
            app.SimLevelLabel_2.HorizontalAlignment = 'right';
            app.SimLevelLabel_2.Layout.Row = 2;
            app.SimLevelLabel_2.Layout.Column = 1;
            app.SimLevelLabel_2.Text = 'Physical Level';

            % Create AltToPlotUnits
            app.AltToPlotUnits = uidropdown(app.GridLayout3);
            app.AltToPlotUnits.Items = {'m', 'hPa'};
            app.AltToPlotUnits.Visible = 'off';
            app.AltToPlotUnits.Layout.Row = 2;
            app.AltToPlotUnits.Layout.Column = 3;
            app.AltToPlotUnits.Value = 'm';

            % Create WarpFlightPathCheckBox
            app.WarpFlightPathCheckBox = uicheckbox(app.GridLayout3);
            app.WarpFlightPathCheckBox.Visible = 'off';
            app.WarpFlightPathCheckBox.Tooltip = {'When choosing different times to plot, warp the flight path to keep it''s location relative to background field.'};
            app.WarpFlightPathCheckBox.Text = 'Warp Flight Path';
            app.WarpFlightPathCheckBox.Layout.Row = 3;
            app.WarpFlightPathCheckBox.Layout.Column = 4;

            % Create WarpTargetsCheckBox
            app.WarpTargetsCheckBox = uicheckbox(app.GridLayout3);
            app.WarpTargetsCheckBox.Visible = 'off';
            app.WarpTargetsCheckBox.Tooltip = {'When choosing different times to plot, warp the targets to keep their location relative to background field.'};
            app.WarpTargetsCheckBox.Text = 'Warp Targets';
            app.WarpTargetsCheckBox.Layout.Row = 4;
            app.WarpTargetsCheckBox.Layout.Column = 4;

            % Create TabGroup
            app.TabGroup = uitabgroup(app.UIFigure);
            app.TabGroup.Position = [1 0 499 517];

            % Create SimulationDetailsTab
            app.SimulationDetailsTab = uitab(app.TabGroup);
            app.SimulationDetailsTab.Title = 'Simulation Details';

            % Create FileSearchPatternEditFieldLabel
            app.FileSearchPatternEditFieldLabel = uilabel(app.SimulationDetailsTab);
            app.FileSearchPatternEditFieldLabel.HorizontalAlignment = 'right';
            app.FileSearchPatternEditFieldLabel.Position = [3 457 107 22];
            app.FileSearchPatternEditFieldLabel.Text = 'File Search Pattern';

            % Create SimulationGlobPattern
            app.SimulationGlobPattern = uieditfield(app.SimulationDetailsTab, 'text');
            app.SimulationGlobPattern.ValueChangedFcn = createCallbackFcn(app, @SimulationGlobPatternValueChanged, true);
            app.SimulationGlobPattern.Tooltip = {''};
            app.SimulationGlobPattern.Placeholder = '/path/to/simulation/directory/*.nc';
            app.SimulationGlobPattern.Position = [123 457 250 22];

            % Create MessagesUI1
            app.MessagesUI1 = uilabel(app.SimulationDetailsTab);
            app.MessagesUI1.VerticalAlignment = 'top';
            app.MessagesUI1.FontSize = 6;
            app.MessagesUI1.Interpreter = 'tex';
            app.MessagesUI1.Position = [14 5 473 210];
            app.MessagesUI1.Text = '';

            % Create FileTimePatternEditFieldLabel
            app.FileTimePatternEditFieldLabel = uilabel(app.SimulationDetailsTab);
            app.FileTimePatternEditFieldLabel.HorizontalAlignment = 'right';
            app.FileTimePatternEditFieldLabel.Position = [14 358 96 22];
            app.FileTimePatternEditFieldLabel.Text = 'File Time Pattern';

            % Create TimeGlobPattern
            app.TimeGlobPattern = uieditfield(app.SimulationDetailsTab, 'text');
            app.TimeGlobPattern.ValueChangedFcn = createCallbackFcn(app, @TimeGlobPatternValueChanged, true);
            app.TimeGlobPattern.Tooltip = {'Enter search pattern for output files. I''ve done my best to guess, correct if necessary. Use the following format: yyMMddHHmmss.'; 'e.g.: ''????????dHHmmss????'''; 'or'; '''?????????sssss????*'''};
            app.TimeGlobPattern.Position = [125 357 248 23];

            % Create FileBrowserButton
            app.FileBrowserButton = uibutton(app.SimulationDetailsTab, 'push');
            app.FileBrowserButton.ButtonPushedFcn = createCallbackFcn(app, @FileBrowserButtonPushed, true);
            app.FileBrowserButton.Position = [383 456 104 23];
            app.FileBrowserButton.Text = 'File Browser';

            % Create SimulationResolutionEditField
            app.SimulationResolutionEditField = uieditfield(app.SimulationDetailsTab, 'numeric');
            app.SimulationResolutionEditField.Limits = [0 Inf];
            app.SimulationResolutionEditField.ValueChangedFcn = createCallbackFcn(app, @SimulationResolutionEditFieldValueChanged, true);
            app.SimulationResolutionEditField.Placeholder = '500';
            app.SimulationResolutionEditField.Position = [125 407 248 22];
            app.SimulationResolutionEditField.Value = 1;

            % Create SimResolutionLabel
            app.SimResolutionLabel = uilabel(app.SimulationDetailsTab);
            app.SimResolutionLabel.HorizontalAlignment = 'right';
            app.SimResolutionLabel.Position = [-10 407 120 22];
            app.SimResolutionLabel.Text = 'Sim. Resolution (m)';

            % Create FlightPlanningTab
            app.FlightPlanningTab = uitab(app.TabGroup);
            app.FlightPlanningTab.Title = 'Flight Planning';

            % Create CalculateBackgroundMotionPanel
            app.CalculateBackgroundMotionPanel = uipanel(app.FlightPlanningTab);
            app.CalculateBackgroundMotionPanel.Title = 'Calculate Background Motion';
            app.CalculateBackgroundMotionPanel.Position = [20 322 189 58];

            % Create CalculateBackgroundMotionButton
            app.CalculateBackgroundMotionButton = uibutton(app.CalculateBackgroundMotionPanel, 'push');
            app.CalculateBackgroundMotionButton.ButtonPushedFcn = createCallbackFcn(app, @CalculateDeformationButtonPushed, true);
            app.CalculateBackgroundMotionButton.Tooltip = {'Calculate the deformation transforms between subsequent time steps using the currently selected variable.'};
            app.CalculateBackgroundMotionButton.Position = [5 9 175 23];
            app.CalculateBackgroundMotionButton.Text = 'Calculate Background Motion';

            % Create ManualFlightPlanningPanel
            app.ManualFlightPlanningPanel = uipanel(app.FlightPlanningTab);
            app.ManualFlightPlanningPanel.Title = 'Manual Flight Planning';
            app.ManualFlightPlanningPanel.Position = [20 390 189 89];

            % Create EditflightlegverticesButton
            app.EditflightlegverticesButton = uibutton(app.ManualFlightPlanningPanel, 'push');
            app.EditflightlegverticesButton.ButtonPushedFcn = createCallbackFcn(app, @EditflightlegverticesButtonPushed, true);
            app.EditflightlegverticesButton.Tooltip = {'Select waypoints in the plot and use plot route to save them to the table. this will also plot the route in the figure.'};
            app.EditflightlegverticesButton.Position = [22 8 144 23];
            app.EditflightlegverticesButton.Text = 'Edit flight leg vertices';

            % Create ManuallyaddflightlegButton
            app.ManuallyaddflightlegButton = uibutton(app.ManualFlightPlanningPanel, 'push');
            app.ManuallyaddflightlegButton.ButtonPushedFcn = createCallbackFcn(app, @ExtractWaypoints, true);
            app.ManuallyaddflightlegButton.Tooltip = {'Select waypoints in the plot and use plot route to save them to the table. this will also plot the route in the figure.'};
            app.ManuallyaddflightlegButton.Position = [22 41 144 23];
            app.ManuallyaddflightlegButton.Text = 'Manually add flight leg';

            % Create AutomaticFlightTargetPlanningPanel
            app.AutomaticFlightTargetPlanningPanel = uipanel(app.FlightPlanningTab);
            app.AutomaticFlightTargetPlanningPanel.Title = 'Automatic Flight/Target Planning';
            app.AutomaticFlightTargetPlanningPanel.Position = [266 357 189 121];

            % Create InteractiveTargetPlannerButton
            app.InteractiveTargetPlannerButton = uibutton(app.AutomaticFlightTargetPlanningPanel, 'push');
            app.InteractiveTargetPlannerButton.ButtonPushedFcn = createCallbackFcn(app, @InteractiveTargetPlannerButtonPushed, true);
            app.InteractiveTargetPlannerButton.Position = [20 71 148 23];
            app.InteractiveTargetPlannerButton.Text = 'Interactive Target Planner';

            % Create InteractiveFlightPlannerButton
            app.InteractiveFlightPlannerButton = uibutton(app.AutomaticFlightTargetPlanningPanel, 'push');
            app.InteractiveFlightPlannerButton.ButtonPushedFcn = createCallbackFcn(app, @InteractiveFlightPlannerButtonPushed, true);
            app.InteractiveFlightPlannerButton.Position = [20 40 148 23];
            app.InteractiveFlightPlannerButton.Text = 'Interactive Flight Planner';

            % Create StartOrContinueFlight
            app.StartOrContinueFlight = uidropdown(app.AutomaticFlightTargetPlanningPanel);
            app.StartOrContinueFlight.Items = {'New', 'Continue', 'Optimize'};
            app.StartOrContinueFlight.Tooltip = {'Start planning a new flight (from scratch), Continue from the end of the existing path, or Optimize the existing path..'};
            app.StartOrContinueFlight.Position = [20 10 148 22];
            app.StartOrContinueFlight.Value = 'New';

            % Create UITable
            app.UITable = uitable(app.FlightPlanningTab);
            app.UITable.ColumnName = {'Start X'; 'Start Y'; 'End X'; 'End Y'; 'Leg Time'; 'Heading'};
            app.UITable.RowName = {};
            app.UITable.Position = [1 115 498 189];

            % Create MessagesUI2
            app.MessagesUI2 = uilabel(app.FlightPlanningTab);
            app.MessagesUI2.VerticalAlignment = 'top';
            app.MessagesUI2.FontSize = 6;
            app.MessagesUI2.Interpreter = 'tex';
            app.MessagesUI2.Position = [14 5 473 111];
            app.MessagesUI2.Text = '';

            % Create AOSPREDetailsTab
            app.AOSPREDetailsTab = uitab(app.TabGroup);
            app.AOSPREDetailsTab.Title = 'AOSPRE Details';

            % Create GridLayout4
            app.GridLayout4 = uigridlayout(app.AOSPREDetailsTab);
            app.GridLayout4.ColumnWidth = {59, 60, 55, 104, 153, '1x'};
            app.GridLayout4.RowHeight = {30, 30, 50, 0, 60, 150, 30, 30, 30};
            app.GridLayout4.ColumnSpacing = 2;
            app.GridLayout4.Padding = [2 10 2 10];

            % Create StartTimeEditFieldLabel
            app.StartTimeEditFieldLabel = uilabel(app.GridLayout4);
            app.StartTimeEditFieldLabel.HorizontalAlignment = 'right';
            app.StartTimeEditFieldLabel.Tooltip = {'select the start time to initialize the simulator.'};
            app.StartTimeEditFieldLabel.Layout.Row = 2;
            app.StartTimeEditFieldLabel.Layout.Column = 2;
            app.StartTimeEditFieldLabel.Text = 'Start Time';

            % Create StartTimeEditField
            app.StartTimeEditField = uieditfield(app.GridLayout4, 'numeric');
            app.StartTimeEditField.Limits = [0 Inf];
            app.StartTimeEditField.ValueDisplayFormat = '%11.0i';
            app.StartTimeEditField.ValueChangedFcn = createCallbackFcn(app, @StartTimeEditFieldValueChanged, true);
            app.StartTimeEditField.Tooltip = {'select the start time to initialize the simulator.'};
            app.StartTimeEditField.Layout.Row = 2;
            app.StartTimeEditField.Layout.Column = [3 4];

            % Create AircraftPanelScanSettings
            app.AircraftPanelScanSettings = uipanel(app.GridLayout4);
            app.AircraftPanelScanSettings.Title = 'Aircraft, Panel, Scan Settings';
            app.AircraftPanelScanSettings.Layout.Row = [5 7];
            app.AircraftPanelScanSettings.Layout.Column = [1 6];

            % Create GridLayout
            app.GridLayout = uigridlayout(app.AircraftPanelScanSettings);
            app.GridLayout.ColumnWidth = {'7.34x', 71, 29, '1x'};
            app.GridLayout.RowHeight = {27, 27, 27, 27, '1x'};
            app.GridLayout.RowSpacing = 9.16666666666667;
            app.GridLayout.Padding = [10 9.16666666666667 10 9.16666666666667];

            % Create Tree
            app.Tree = uitree(app.GridLayout, 'checkbox');
            app.Tree.Tooltip = {'Select which panels/scans to use. '};
            app.Tree.Layout.Row = [1 5];
            app.Tree.Layout.Column = 1;

            % Create LHSPortPanelNode
            app.LHSPortPanelNode = uitreenode(app.Tree);
            app.LHSPortPanelNode.Text = 'LHS/Port Panel';

            % Create LHSRHI
            app.LHSRHI = uitreenode(app.LHSPortPanelNode);
            app.LHSRHI.Text = 'RHI';

            % Create LHSPPI
            app.LHSPPI = uitreenode(app.LHSPortPanelNode);
            app.LHSPPI.Text = 'PPI';

            % Create LHSFore
            app.LHSFore = uitreenode(app.LHSPortPanelNode);
            app.LHSFore.Text = 'Fore';

            % Create AftNode_2
            app.AftNode_2 = uitreenode(app.LHSPortPanelNode);
            app.AftNode_2.Text = 'Aft';

            % Create RHSStarboardPanelNode
            app.RHSStarboardPanelNode = uitreenode(app.Tree);
            app.RHSStarboardPanelNode.Text = 'RHS/Starboard Panel';

            % Create RHSRHI
            app.RHSRHI = uitreenode(app.RHSStarboardPanelNode);
            app.RHSRHI.Text = 'RHI';

            % Create RHSPPI
            app.RHSPPI = uitreenode(app.RHSStarboardPanelNode);
            app.RHSPPI.Text = 'PPI';

            % Create RHSFore
            app.RHSFore = uitreenode(app.RHSStarboardPanelNode);
            app.RHSFore.Text = 'Fore';

            % Create AftNode
            app.AftNode = uitreenode(app.RHSStarboardPanelNode);
            app.AftNode.Text = 'Aft';

            % Create TopPanelNode
            app.TopPanelNode = uitreenode(app.Tree);
            app.TopPanelNode.Text = 'Top Panel';

            % Create TOPPPI
            app.TOPPPI = uitreenode(app.TopPanelNode);
            app.TOPPPI.Text = 'Aft';

            % Create TOPFore
            app.TOPFore = uitreenode(app.TopPanelNode);
            app.TOPFore.Text = 'Fore';

            % Create BottomPanelNode
            app.BottomPanelNode = uitreenode(app.Tree);
            app.BottomPanelNode.Text = 'Bottom Panel';

            % Create BOTFore
            app.BOTFore = uitreenode(app.BottomPanelNode);
            app.BOTFore.Text = 'Fore';

            % Create AftNode_3
            app.AftNode_3 = uitreenode(app.BottomPanelNode);
            app.AftNode_3.Text = 'Aft';

            % Assign Checked Nodes
            app.Tree.CheckedNodesChangedFcn = createCallbackFcn(app, @GenerateNamelists, true);

            % Create AngularResolutionLabel
            app.AngularResolutionLabel = uilabel(app.GridLayout);
            app.AngularResolutionLabel.HorizontalAlignment = 'right';
            app.AngularResolutionLabel.Layout.Row = 2;
            app.AngularResolutionLabel.Layout.Column = [2 3];
            app.AngularResolutionLabel.Text = 'Angular Resolution';

            % Create AngularResolutionEditField
            app.AngularResolutionEditField = uieditfield(app.GridLayout, 'numeric');
            app.AngularResolutionEditField.ValueChangedFcn = createCallbackFcn(app, @RecalculateScanlistsValueChanged, true);
            app.AngularResolutionEditField.Layout.Row = 2;
            app.AngularResolutionEditField.Layout.Column = 4;
            app.AngularResolutionEditField.Value = 1.5;

            % Create MaximumTiltRotLabel
            app.MaximumTiltRotLabel = uilabel(app.GridLayout);
            app.MaximumTiltRotLabel.HorizontalAlignment = 'right';
            app.MaximumTiltRotLabel.Visible = 'off';
            app.MaximumTiltRotLabel.Layout.Row = 5;
            app.MaximumTiltRotLabel.Layout.Column = [2 3];
            app.MaximumTiltRotLabel.Text = 'Maximum Tilt/Rot';

            % Create MaximumTiltRotEditField
            app.MaximumTiltRotEditField = uieditfield(app.GridLayout, 'numeric');
            app.MaximumTiltRotEditField.ValueChangedFcn = createCallbackFcn(app, @RecalculateScanlistsValueChanged, true);
            app.MaximumTiltRotEditField.Visible = 'off';
            app.MaximumTiltRotEditField.Layout.Row = 5;
            app.MaximumTiltRotEditField.Layout.Column = 4;
            app.MaximumTiltRotEditField.Value = 45;

            % Create ForeEditFieldLabel
            app.ForeEditFieldLabel = uilabel(app.GridLayout);
            app.ForeEditFieldLabel.HorizontalAlignment = 'right';
            app.ForeEditFieldLabel.Layout.Row = 3;
            app.ForeEditFieldLabel.Layout.Column = [2 3];
            app.ForeEditFieldLabel.Text = 'Fore';

            % Create ForeEditField
            app.ForeEditField = uieditfield(app.GridLayout, 'numeric');
            app.ForeEditField.ValueChangedFcn = createCallbackFcn(app, @RecalculateScanlistsValueChanged, true);
            app.ForeEditField.Layout.Row = 3;
            app.ForeEditField.Layout.Column = 4;
            app.ForeEditField.Value = 5;

            % Create AftEditFieldLabel
            app.AftEditFieldLabel = uilabel(app.GridLayout);
            app.AftEditFieldLabel.HorizontalAlignment = 'right';
            app.AftEditFieldLabel.Layout.Row = 4;
            app.AftEditFieldLabel.Layout.Column = [2 3];
            app.AftEditFieldLabel.Text = 'Aft';

            % Create AftEditField
            app.AftEditField = uieditfield(app.GridLayout, 'numeric');
            app.AftEditField.ValueChangedFcn = createCallbackFcn(app, @RecalculateScanlistsValueChanged, true);
            app.AftEditField.Layout.Row = 4;
            app.AftEditField.Layout.Column = 4;
            app.AftEditField.Value = 35;

            % Create AircraftVelmsEditFieldLabel
            app.AircraftVelmsEditFieldLabel = uilabel(app.GridLayout);
            app.AircraftVelmsEditFieldLabel.HorizontalAlignment = 'right';
            app.AircraftVelmsEditFieldLabel.Layout.Row = 1;
            app.AircraftVelmsEditFieldLabel.Layout.Column = [2 3];
            app.AircraftVelmsEditFieldLabel.Text = 'Aircraft Vel. (m/s)';

            % Create AircraftVelocity
            app.AircraftVelocity = uieditfield(app.GridLayout, 'numeric');
            app.AircraftVelocity.Layout.Row = 1;
            app.AircraftVelocity.Layout.Column = 4;
            app.AircraftVelocity.Value = 120;

            % Create pathtoAOSPRELabel
            app.pathtoAOSPRELabel = uilabel(app.GridLayout4);
            app.pathtoAOSPRELabel.HorizontalAlignment = 'right';
            app.pathtoAOSPRELabel.Layout.Row = 1;
            app.pathtoAOSPRELabel.Layout.Column = [1 2];
            app.pathtoAOSPRELabel.Text = 'path to AOSPRE';

            % Create AOSPREPath
            app.AOSPREPath = uieditfield(app.GridLayout4, 'text');
            app.AOSPREPath.Tooltip = {''};
            app.AOSPREPath.Placeholder = '/path/to/AOSPRE';
            app.AOSPREPath.Layout.Row = 1;
            app.AOSPREPath.Layout.Column = [3 5];

            % Create StartTimeSlider
            app.StartTimeSlider = uislider(app.GridLayout4);
            app.StartTimeSlider.ValueChangedFcn = createCallbackFcn(app, @StartTimeSliderValueChanged, true);
            app.StartTimeSlider.FontSize = 8;
            app.StartTimeSlider.Layout.Row = 3;
            app.StartTimeSlider.Layout.Column = [1 6];

            % Create SameasplottedtimeCheckBox
            app.SameasplottedtimeCheckBox = uicheckbox(app.GridLayout4);
            app.SameasplottedtimeCheckBox.Text = 'Same as plotted time';
            app.SameasplottedtimeCheckBox.Layout.Row = 2;
            app.SameasplottedtimeCheckBox.Layout.Column = [5 6];
            app.SameasplottedtimeCheckBox.Value = true;

            % Create InitAOSPREButton
            app.InitAOSPREButton = uibutton(app.GridLayout4, 'push');
            app.InitAOSPREButton.ButtonPushedFcn = createCallbackFcn(app, @InitAOSPREButtonPushed, true);
            app.InitAOSPREButton.Layout.Row = 8;
            app.InitAOSPREButton.Layout.Column = [5 6];
            app.InitAOSPREButton.Text = 'Init. AOSPRE';

            % Create MessagesUI3
            app.MessagesUI3 = uilabel(app.GridLayout4);
            app.MessagesUI3.VerticalAlignment = 'top';
            app.MessagesUI3.FontSize = 6;
            app.MessagesUI3.Layout.Row = [8 9];
            app.MessagesUI3.Layout.Column = [1 4];
            app.MessagesUI3.Interpreter = 'tex';
            app.MessagesUI3.Text = '';

            % Create AdvancedOptionsTab
            app.AdvancedOptionsTab = uitab(app.TabGroup);
            app.AdvancedOptionsTab.Title = 'Advanced Options';

            % Create outputformatstringLabel
            app.outputformatstringLabel = uilabel(app.AdvancedOptionsTab);
            app.outputformatstringLabel.HorizontalAlignment = 'right';
            app.outputformatstringLabel.WordWrap = 'on';
            app.outputformatstringLabel.Position = [4 456 176 32];
            app.outputformatstringLabel.Text = 'output_filename_format_string = ';

            % Create output_filename_format_string
            app.output_filename_format_string = uieditfield(app.AdvancedOptionsTab, 'text');
            app.output_filename_format_string.Tooltip = {'The naming pattern to use on the outputted CFradial netCDF files. "A" is the placeholder to AOSPRE to input time details.'; ''; 'PNL (SCN) will be replaced with the panel (scan) identifier.'};
            app.output_filename_format_string.Position = [149 437 246 22];
            app.output_filename_format_string.Value = '''("./output/PNL_SCN_",A,"_to_",A,".nc")''';

            % Create flight_level_coordinateEditFieldLabel
            app.flight_level_coordinateEditFieldLabel = uilabel(app.AdvancedOptionsTab);
            app.flight_level_coordinateEditFieldLabel.HorizontalAlignment = 'right';
            app.flight_level_coordinateEditFieldLabel.WordWrap = 'on';
            app.flight_level_coordinateEditFieldLabel.Position = [5 398 176 32];
            app.flight_level_coordinateEditFieldLabel.Text = 'flight_level_coordinate = ';

            % Create flight_level_coordinate
            app.flight_level_coordinate = uieditfield(app.AdvancedOptionsTab, 'text');
            app.flight_level_coordinate.Tooltip = {'The naming pattern to use on the outputted CFradial netCDF files. "A" is the placeholder to AOSPRE to input time details.'; ''; 'PNL (SCN) will be replaced with the panel (scan) identifier.'};
            app.flight_level_coordinate.Position = [150 379 246 22];
            app.flight_level_coordinate.Value = '"Z"';

            % Create flight_level_waypoints_vertEditFieldLabel
            app.flight_level_waypoints_vertEditFieldLabel = uilabel(app.AdvancedOptionsTab);
            app.flight_level_waypoints_vertEditFieldLabel.HorizontalAlignment = 'right';
            app.flight_level_waypoints_vertEditFieldLabel.WordWrap = 'on';
            app.flight_level_waypoints_vertEditFieldLabel.Position = [4 341 176 32];
            app.flight_level_waypoints_vertEditFieldLabel.Text = 'flight_level_waypoints_vert = ';

            % Create flight_level_waypoints_vert
            app.flight_level_waypoints_vert = uieditfield(app.AdvancedOptionsTab, 'text');
            app.flight_level_waypoints_vert.Tooltip = {'The naming pattern to use on the outputted CFradial netCDF files. "A" is the placeholder to AOSPRE to input time details.'; ''; 'PNL (SCN) will be replaced with the panel (scan) identifier.'};
            app.flight_level_waypoints_vert.Position = [149 322 246 22];
            app.flight_level_waypoints_vert.Value = '2000';

            % Create ScanningOptionsPanel
            app.ScanningOptionsPanel = uipanel(app.AdvancedOptionsTab);
            app.ScanningOptionsPanel.Tooltip = {'Some advanced '};
            app.ScanningOptionsPanel.Title = 'Scanning Options';
            app.ScanningOptionsPanel.Position = [19 5 364 297];

            % Create ScanningOptionsGrid
            app.ScanningOptionsGrid = uigridlayout(app.ScanningOptionsPanel);
            app.ScanningOptionsGrid.ColumnWidth = {39, 98, 36, '1x'};
            app.ScanningOptionsGrid.RowHeight = {23, 22, 22, 22, 22, 22, 22, 22, 22};
            app.ScanningOptionsGrid.ColumnSpacing = 9.4;
            app.ScanningOptionsGrid.RowSpacing = 7.83333333333333;
            app.ScanningOptionsGrid.Padding = [9.4 7.83333333333333 9.4 7.83333333333333];

            % Create pulse_repetition_frequencyEditFieldLabel
            app.pulse_repetition_frequencyEditFieldLabel = uilabel(app.ScanningOptionsGrid);
            app.pulse_repetition_frequencyEditFieldLabel.HorizontalAlignment = 'right';
            app.pulse_repetition_frequencyEditFieldLabel.Layout.Row = 2;
            app.pulse_repetition_frequencyEditFieldLabel.Layout.Column = [1 3];
            app.pulse_repetition_frequencyEditFieldLabel.Text = 'pulse_repetition_frequency =';

            % Create pulse_repetition_frequency
            app.pulse_repetition_frequency = uieditfield(app.ScanningOptionsGrid, 'numeric');
            app.pulse_repetition_frequency.Layout.Row = 2;
            app.pulse_repetition_frequency.Layout.Column = 4;
            app.pulse_repetition_frequency.Value = 2500;

            % Create pulses_per_pulse_setEditFieldLabel
            app.pulses_per_pulse_setEditFieldLabel = uilabel(app.ScanningOptionsGrid);
            app.pulses_per_pulse_setEditFieldLabel.HorizontalAlignment = 'right';
            app.pulses_per_pulse_setEditFieldLabel.Layout.Row = 3;
            app.pulses_per_pulse_setEditFieldLabel.Layout.Column = [2 3];
            app.pulses_per_pulse_setEditFieldLabel.Text = 'pulses_per_pulse_set = ';

            % Create pulses_per_pulse_set
            app.pulses_per_pulse_set = uieditfield(app.ScanningOptionsGrid, 'numeric');
            app.pulses_per_pulse_set.Layout.Row = 3;
            app.pulses_per_pulse_set.Layout.Column = 4;
            app.pulses_per_pulse_set.Value = 2;

            % Create revisits_per_acquisition_timeEditFieldLabel
            app.revisits_per_acquisition_timeEditFieldLabel = uilabel(app.ScanningOptionsGrid);
            app.revisits_per_acquisition_timeEditFieldLabel.HorizontalAlignment = 'right';
            app.revisits_per_acquisition_timeEditFieldLabel.Layout.Row = 4;
            app.revisits_per_acquisition_timeEditFieldLabel.Layout.Column = [1 3];
            app.revisits_per_acquisition_timeEditFieldLabel.Text = 'revisits_per_acquisition_time = ';

            % Create revisits_per_acquisition_time
            app.revisits_per_acquisition_time = uieditfield(app.ScanningOptionsGrid, 'numeric');
            app.revisits_per_acquisition_time.Layout.Row = 4;
            app.revisits_per_acquisition_time.Layout.Column = 4;
            app.revisits_per_acquisition_time.Value = 20;

            % Create beams_per_acquisition_timeEditFieldLabel
            app.beams_per_acquisition_timeEditFieldLabel = uilabel(app.ScanningOptionsGrid);
            app.beams_per_acquisition_timeEditFieldLabel.HorizontalAlignment = 'right';
            app.beams_per_acquisition_timeEditFieldLabel.Layout.Row = 5;
            app.beams_per_acquisition_timeEditFieldLabel.Layout.Column = [1 3];
            app.beams_per_acquisition_timeEditFieldLabel.Text = 'beams_per_acquisition_time';

            % Create beams_per_acquisition_time
            app.beams_per_acquisition_time = uieditfield(app.ScanningOptionsGrid, 'numeric');
            app.beams_per_acquisition_time.Layout.Row = 5;
            app.beams_per_acquisition_time.Layout.Column = 4;
            app.beams_per_acquisition_time.Value = 6;

            % Create CRSIM_ConfigPath
            app.CRSIM_ConfigPath = uieditfield(app.ScanningOptionsGrid, 'text');
            app.CRSIM_ConfigPath.Layout.Row = 1;
            app.CRSIM_ConfigPath.Layout.Column = 4;
            app.CRSIM_ConfigPath.Value = '"CONFIG_crsim"';

            % Create CRSIM_ConfigEditFieldLabel
            app.CRSIM_ConfigEditFieldLabel = uilabel(app.ScanningOptionsGrid);
            app.CRSIM_ConfigEditFieldLabel.HorizontalAlignment = 'right';
            app.CRSIM_ConfigEditFieldLabel.Layout.Row = 1;
            app.CRSIM_ConfigEditFieldLabel.Layout.Column = [2 3];
            app.CRSIM_ConfigEditFieldLabel.Text = 'CRSIM_Config = ';

            % Create meters_between_gatesEditFieldLabel
            app.meters_between_gatesEditFieldLabel = uilabel(app.ScanningOptionsGrid);
            app.meters_between_gatesEditFieldLabel.HorizontalAlignment = 'right';
            app.meters_between_gatesEditFieldLabel.Layout.Row = 6;
            app.meters_between_gatesEditFieldLabel.Layout.Column = [2 3];
            app.meters_between_gatesEditFieldLabel.Text = 'meters _between_gates = ';

            % Create meters_between_gates
            app.meters_between_gates = uieditfield(app.ScanningOptionsGrid, 'text');
            app.meters_between_gates.HorizontalAlignment = 'right';
            app.meters_between_gates.Layout.Row = 6;
            app.meters_between_gates.Layout.Column = 4;
            app.meters_between_gates.Value = '75';

            % Create meters_to_center_of_first_gateEditFieldLabel
            app.meters_to_center_of_first_gateEditFieldLabel = uilabel(app.ScanningOptionsGrid);
            app.meters_to_center_of_first_gateEditFieldLabel.HorizontalAlignment = 'right';
            app.meters_to_center_of_first_gateEditFieldLabel.Layout.Row = 7;
            app.meters_to_center_of_first_gateEditFieldLabel.Layout.Column = [1 3];
            app.meters_to_center_of_first_gateEditFieldLabel.Text = 'meters_to_center_of_first_gate = ';

            % Create meters_to_center_of_first_gate
            app.meters_to_center_of_first_gate = uieditfield(app.ScanningOptionsGrid, 'numeric');
            app.meters_to_center_of_first_gate.Layout.Row = 7;
            app.meters_to_center_of_first_gate.Layout.Column = 4;
            app.meters_to_center_of_first_gate.Value = 1;

            % Create max_range_in_metersEditFieldLabel
            app.max_range_in_metersEditFieldLabel = uilabel(app.ScanningOptionsGrid);
            app.max_range_in_metersEditFieldLabel.HorizontalAlignment = 'right';
            app.max_range_in_metersEditFieldLabel.Layout.Row = 8;
            app.max_range_in_metersEditFieldLabel.Layout.Column = [2 3];
            app.max_range_in_metersEditFieldLabel.Text = 'max_range_in_meters = ';

            % Create max_range_in_meters
            app.max_range_in_meters = uieditfield(app.ScanningOptionsGrid, 'numeric');
            app.max_range_in_meters.Layout.Row = 8;
            app.max_range_in_meters.Layout.Column = 4;
            app.max_range_in_meters.Value = 75000;

            % Create skip_seconds_between_scansLabel
            app.skip_seconds_between_scansLabel = uilabel(app.ScanningOptionsGrid);
            app.skip_seconds_between_scansLabel.HorizontalAlignment = 'right';
            app.skip_seconds_between_scansLabel.Layout.Row = 9;
            app.skip_seconds_between_scansLabel.Layout.Column = [2 3];
            app.skip_seconds_between_scansLabel.Text = 'skip_seconds_between_scans = ';

            % Create skip_seconds_between_scans
            app.skip_seconds_between_scans = uieditfield(app.ScanningOptionsGrid, 'numeric');
            app.skip_seconds_between_scans.Layout.Row = 9;
            app.skip_seconds_between_scans.Layout.Column = 4;

            % Create ContextMenu
            app.ContextMenu = uicontextmenu(app.UIFigure);

            % Create ClearMenu
            app.ClearMenu = uimenu(app.ContextMenu);
            app.ClearMenu.Text = 'Clear';

            % Create ClearSelection
            app.ClearSelection = uimenu(app.ClearMenu);
            app.ClearSelection.MenuSelectedFcn = createCallbackFcn(app, @ClearSelectionSelected, true);
            app.ClearSelection.Text = 'Selection';

            % Create LegRowMenu
            app.LegRowMenu = uimenu(app.ClearMenu);
            app.LegRowMenu.MenuSelectedFcn = createCallbackFcn(app, @LegRowMenuSelected, true);
            app.LegRowMenu.Text = 'Leg/Row';

            % Create RefreshMenu
            app.RefreshMenu = uimenu(app.ContextMenu);
            app.RefreshMenu.Text = 'Refresh';

            % Create RefreshSelectionMenu
            app.RefreshSelectionMenu = uimenu(app.RefreshMenu);
            app.RefreshSelectionMenu.MenuSelectedFcn = createCallbackFcn(app, @RefreshSelectionMenuSelected, true);
            app.RefreshSelectionMenu.Text = 'Selection';

            % Create AddnewLegMenu
            app.AddnewLegMenu = uimenu(app.ContextMenu);
            app.AddnewLegMenu.MenuSelectedFcn = createCallbackFcn(app, @AddnewLegMenuSelected, true);
            app.AddnewLegMenu.Text = 'Add new Leg';
            
            % Assign app.ContextMenu
            app.UITable.ContextMenu = app.ContextMenu;

            % Create ContextMenu2
            app.ContextMenu2 = uicontextmenu(app.UIFigure);

            % Create RefreshMenu_2
            app.RefreshMenu_2 = uimenu(app.ContextMenu2);
            app.RefreshMenu_2.MenuSelectedFcn = createCallbackFcn(app, @RefreshMenu_2Selected, true);
            app.RefreshMenu_2.Text = 'Refresh';
            
            % Assign app.ContextMenu2
            app.UIAxes.ContextMenu = app.ContextMenu2;

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = aospregui

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFunction)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end