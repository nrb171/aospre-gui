classdef aospregui < matlab.apps.AppBase

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
        FileSearchFeedback              matlab.ui.control.Label
        SimulationGlobPattern           matlab.ui.control.EditField
        FileSearchPatternEditFieldLabel  matlab.ui.control.Label
        FlightPlanningTab               matlab.ui.container.Tab
        EditpathButton                  matlab.ui.control.Button
        StartOrContinueFlight           matlab.ui.control.DropDown
        InteractiveFlightPlannerButton  matlab.ui.control.Button
        UITable                         matlab.ui.control.Table
        ManuallydrawpathButton          matlab.ui.control.Button
        InteractiveTargetPlannerButton  matlab.ui.control.Button
        AOSPREDetailsTab                matlab.ui.container.Tab
        GridLayout4                     matlab.ui.container.GridLayout
        StartTimeSlider                 matlab.ui.control.Slider
        AOSPREPath                      matlab.ui.control.EditField
        pathtoAOSPRELabel               matlab.ui.control.Label
        InitAOSPREButton                matlab.ui.control.Button
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
        max_range_in_meters             matlab.ui.control.NumericEditField
        max_range_in_metersEditFieldLabel  matlab.ui.control.Label
        meters_to_center_of_first_gate  matlab.ui.control.NumericEditField
        meters_to_center_of_first_gateEditFieldLabel  matlab.ui.control.Label
        meters_between_gatesEditField   matlab.ui.control.EditField
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
        PlottedVariableDropDown         matlab.ui.control.DropDown
        PlottedVariableDropDownLabel    matlab.ui.control.Label
        SimLevelLabel                   matlab.ui.control.Label
        SimulationLevelEditField        matlab.ui.control.NumericEditField
        SimulationLevelSlider           matlab.ui.control.Slider
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
        
        environment           % array of environment vertical slice.
    end
    properties (Access = public)
        flightPath          % roi polyline of flight path.     
        Targets = struct(...
            'map',[], ...    % assigned after region finding. 2D map. Values: 0=no target, integer = targetID.
            'weight',[], ... % default 10, can be set
            'dwellTime', []) % default 0 (no limit), can be set.

    end
    
    methods (Access = private)
        
        function [] = UpdateFigure(app)
            % function to redraw plan view

            
            %% get data
            var = ncread(...
                [...            % load variable at appropriate time
                    app.FileList(app.TimeIndexEditField.Value).folder, '/', ...
                    app.FileList(app.TimeIndexEditField.Value).name ...
                ], ...
                app.PlottedVariableDropDown.Value ...   %variable to load
            );
            
            varSlice = var(:,:,app.SimulationLevelEditField.Value);
            app.environment = varSlice';

            %% draw map

            xlat = ncread(...
                [...            % load variable at appropriate time
                    app.FileList(app.TimeIndexEditField.Value).folder, '/', ...
                    app.FileList(app.TimeIndexEditField.Value).name ...
                ], ...
                'XLAT' ...   %variable to load
            );
            xlon = ncread(...
                [...            % load variable at appropriate time
                    app.FileList(app.TimeIndexEditField.Value).folder, '/', ...
                    app.FileList(app.TimeIndexEditField.Value).name ...
                ], ...
                'XLONG' ...   %variable to load
            );






            imagesc(app.UIAxes, ...
                [min(xlon(:)), max(xlon(:))], ...
                [min(xlat(:)), max(xlat(:))], ...
                var ...
            )


            p = pcolor(app.UIAxes, varSlice'); % transpose to fix WRF/MATLAB row/column difference

            p.EdgeColor = 'None';
            colorbar(app.UIAxes)
            clim(app.UIAxes, [prctile(varSlice(:), 1), prctile(varSlice(:), 99)+1e-3])

            %% draw path if one exists
            if ~isempty(app.UITable.Data)
                hold(app.UIAxes, 'on')

                x0 = app.UITable.Data(:,[1,3]);
                x = [x0(1,1), x0(:,2)'];
                y0 = app.UITable.Data(:,[2,4]);
                y = [y0(1,1), y0(:,2)'];

                for i = 1:numel(x)-1
                    p1 = [x(i) y(i)];                         % First Point
                    p2 = [x(i+1) y(i+1)];                         % Second Point
                    dp = p2-p1;                         % Difference

                    h=quiver(app.UIAxes,p1(1),p1(2),dp(1),dp(2),0, 'filled', 'r', 'LineWidth',1);
                    
                    %set(h,'MaxHeadSize',1,'AutoScaleFactor',1);
                    
                end
                scatter(app.UIAxes, x, y, 'r', 'filled')
                
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
                    contour(app.UIAxes, app.Targets.map == uTargets(i), 1, 'LineColor', cmap(i,:), 'LineWidth', 2)
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
        
        function [] = ChildrenLevel2Visibility(app)
            app.PlottedVariableDropDownLabel.Visible = "on";
            app.SimulationLevelSlider.Visible = "on";
            app.SimLevelLabel.Visible = "on";
            app.SimulationLevelEditField.Visible = "on";
            app.TimeIndexEditField.Visible = "on";
            app.TimeIndexLabel.Visible = "on";
            app.TimeIndexSlider.Visible = "on";
            app.PlottedVariableDropDown.Visible = "on";
            app.SimResolutionLabel.Visible = "on";
            app.SimulationResolutionEditField.Visible = "on";
            
        end
        
        
        function [] = UpdateSlidersAndRanges(app)
            %% set simulation level 
            % get temporary variable to set number of levels/times
            tempVar = ncread(...
                [...    %filename
                    app.FileList(1).folder,'/', ...
                    app.FileList(1).name...
                ], ...  %variable name
                app.PlottedVariableDropDown.Value ...
            );
            
            
            % get levels/times
            verticalLevels = size(tempVar,3);
            timeIndices = numel(app.FileList);
            
            % set ranges and temp values for sliders/edit fields
            app.TimeIndexEditField.Value = 1;
            app.TimeIndexSlider.Limits = [1, timeIndices];
            app.TimeIndexSlider.MinorTicks = 1:ceil(log10(timeIndices)):timeIndices;
            app.TimeIndexSlider.MajorTicks = intersect(round(linspace(1, timeIndices, min(sqrt(timeIndices), 10))), app.TimeIndexSlider.MinorTicks);
            app.TimeIndexSlider.Value = 1;
            

            [~,bestLevelGuess] = max(squeeze(std(std(tempVar))));
            app.SimulationLevelEditField.Value = bestLevelGuess;
            app.SimulationLevelSlider.Value = bestLevelGuess;
            app.SimulationLevelSlider.Limits = [1, verticalLevels];
            app.SimulationLevelSlider.MinorTicks = 1:verticalLevels;
            app.SimulationLevelSlider.MajorTicks = intersect(round(linspace(1, verticalLevels, sqrt(verticalLevels))), app.SimulationLevelSlider.MinorTicks);
        
        end
        
        function [] = UpdateTableLogic(app)
            % try to update the tables based on some simple logic

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
        
        
        
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            %% add helpers folder to path
            addpath(genpath("./helpers/"))
            
            %% make .meta/* or output folder if none exists
            if ~exist("./.meta", "file")
                mkdir("./.meta")
            end

            if ~exist("./.meta/scanTables", "file")
                mkdir("./.meta/scanTables")
            end

            if ~exist("./output", "file")
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
            if exist("./.meta/aospregui.log","file")
                FID = fopen('./.meta/aospregui.log');
                app.SimulationGlobPattern.Value = fgetl(FID);
                app.TimeGlobPattern.Value = fgetl(FID);
                app.SimulationResolutionEditField.Value = str2double(fgetl(FID));
                app.AOSPREPath.Value = fgetl(FID);
            end
            if exist("./.meta/table.txt", 'file')
                TAB = readmatrix('./.meta/table.txt');
                app.UITable.Data = TAB;
            end
            app.SimulationGlobPatternValueChanged(app)
            app.TimeGlobPatternValueChanged(app)
        end

        % Value changed function: StartTimeEditField
        function StartTimeEditFieldValueChanged(app, event)
            value = app.StartTimeEditField.Value;
            [~,nearestIdx] = min(abs(value - app.StartTimeSlider.MajorTicks));
            value = app.StartTimeSlider.MajorTicks(nearestIdx);
            app.StartTimeEditField.Value = value;
            app.StartTimeSlider.Value = value;
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
                    app.FileSearchFeedback.Text = num2str(numel(app.FileList))+" files found.";
    
                
                    %% update variable dropdown based on file search
                    % get variables in file
                    
                    simInfo=ncinfo([app.FileList(1).folder,'/',app.FileList(1).name],"/");
                    varNames = string({simInfo.Variables.Name});
                    app.PlottedVariableDropDown.Items = varNames;
                    
                    %% Set children visibility
                    ChildrenLevel2Visibility(app)

                    

                    %% Try variable    
                    try
                        
                        varsToCheck = ["W", "Q", "QVAPOR"];
                        for iVar = 1:numel(varsToCheck)
                            
                            mask = strcmp(varNames, varsToCheck(iVar));
                            if nnz(mask) > 0
                                app.PlottedVariableDropDown.Value = varNames(find(mask, 1));
                                app.UpdateSlidersAndRanges()
                                app.PlottedVariableDropDownValueChanged()
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
                       commonResolutions = [50, 100, 250, 500, 1000, 3000, 9000, 27000]; % m
                       dlat = XLAT(2,2,1)-XLAT(1,1,1);
                       dy = dlat*111131.745;
                       
                       [~,bestGuessResolutionIdx] = min(abs(dy-commonResolutions));
                       bestGuessResolution = commonResolutions(bestGuessResolutionIdx);
                       app.SimulationResolutionEditField.Value = bestGuessResolution;
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
                        app.TabGroup.SelectedTab = app.FlightPlanningTab;
                    end

                else
                    app.FileSearchFeedback.Text = "Choose directory with at least 2 files.";
                end
            catch ME
                ME.message
                app.FileSearchFeedback.Text = {'Search Failed, try again.', char(ME.message)};
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

        % Value changed function: PlottedVariableDropDown
        function PlottedVariableDropDownValueChanged(app, event)
            value = app.PlottedVariableDropDown.Value;
            % update sliders
            UpdateSlidersAndRanges(app)
            % update figure
            UpdateFigure(app);
       
        end

        % Value changed function: TimeIndexSlider
        function TimeIndexSliderValueChanged(app, event)
            % link edit field to slider value
            value = app.TimeIndexSlider.Value;
            app.TimeIndexSlider.Value = round(value);
            app.TimeIndexEditField.Value = round(value);

            % update figure
            UpdateFigure(app);
        end

        % Value changed function: TimeIndexEditField
        function TimeIndexEditFieldValueChanged(app, event)
            % link slider value to edit field
            value = app.TimeIndexEditField.Value;
            if value >= min(app.TimeIndexSlider.Limits) && ...
                    value <= max(app.TimeIndexSlider.Limits)
                app.TimeIndexSlider.Value = round(value);

                % update figure
                UpdateFigure(app);
            else
                app.TimeIndexEditField.Value = 1;
            end
        end

        % Value changed function: SimulationLevelSlider
        function SimulationLevelSliderValueChanged(app, event)
            %link edit field to slider value
            value = app.SimulationLevelSlider.Value;
            app.SimulationLevelSlider.Value = round(value);
            app.SimulationLevelEditField.Value = round(value);

            % update figure
            UpdateFigure(app);
        end

        % Value changed function: SimulationLevelEditField
        function SimulationLevelEditFieldValueChanged(app, event)
            value = app.SimulationLevelEditField.Value;
            if value >= min(app.SimulationLevelSlider.Limits) && ...
                    value <= max(app.SimulationLevelSlider.Limits)
                app.SimulationLevelSlider.Value = round(value);

                % update figure
                UpdateFigure(app);
            else
                app.SimulationLevelEditField.Value = 1;
            end
                
        end

        % Button pushed function: ManuallydrawpathButton
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
                yc = y(~isnan(x) | ~isnan(y));
                app.flightPath = images.roi.Polyline(app.UIAxes);
                beginDrawingFromPoint(app.flightPath, [xc(end), yc(end)])
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

        % Cell edit callback: UITable
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
            %% link StartTimeEditField
            [~, ind] = min(abs(value(1) - app.FileTimes));
            app.StartTimeSlider.Value = app.FileTimes(ind);
            app.StartTimeEditField.Value = app.FileTimes(ind);
            
            
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

            value = app.AngularResolutionEditField.Value;

            % call the functions
            GenerateScanlists(app)
            SaveScanlists(app)

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
            checkedNodes = app.Tree.CheckedNodes;

            % clear all prio namelists
            delete("./output/namelists/*.txt")
            delete("./output/scanlists/*.txt")

            %% loop through all nodes, finding the parent nodes.
            files = [];
            for iNode = 1:numel(checkedNodes)
                if numel(checkedNodes(iNode).Parent)>0
                    try
                        if ~exist("./output/namelists/"+lower(checkedNodes(iNode).Parent.Text(1:3)+"_namelist.txt"), 'file')
                            WriteNamelist(app, checkedNodes(iNode))
                        end
                        
                        % add relevant scan table to final scanlist.
                        fileID1 = fopen(sprintf("./.meta/scanTables/%s_%s.txt",lower(checkedNodes(iNode).Parent.Text(1:3)), lower(checkedNodes(iNode).Text)), 'r');
                        fname = sprintf("./output/scanlists/%s.txt",lower(checkedNodes(iNode).Parent.Text(1:3)));
                    
                        
                         
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
                        fclose(fileID2)
                        fclose(fileID1)
                    catch ME
                        if isa(ME.identifier, 'MATLAB:FileIO:InvalidFid')
                            continue
                        end
                    end
                       
                end
            end

            function WriteNamelist(app, node)
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
                try
                    namelistText(end + 1) = " leg_initial_time = " + app.StartTimeEditField.Value;
                    namelistText(end + 1) = " leg_time_seconds = " + strjoin(string(floor(TAB(:,5))), ", ");
                    
                    namelistText(end + 1) = " time_evolution = .TRUE.";
                    
                    namelistText(end + 1) = " flight_waypoints_x = " + ...
                        strjoin(...
                            string(...
                                TABT([1,3:6:3+6*(size(TAB,1)-1)])...
                            ), ...
                        ", "...
                    );
                    namelistText(end + 1) = " flight_waypoints_y = " + ...
                        strjoin(...
                            string(...
                                TABT([2,4:6:4+6*(size(TAB,1)-1)])...
                            ), ...
                        ", " ...
                    );
                    namelistText(end + 1) = " flight_waypoints_vert = " + ...
                        strjoin(...
                            string(...
                                ones(size(TABT([2,4:6:4+6*(size(TAB,1)-1)])))*str2double(app.flight_level_waypoints_vert.Value)...
                            ), ...
                        ", " ...
                    );

                    namelistText(end + 1) =  " flight_level_coordinate = "+app.flight_level_coordinate.Value;
                catch
                    error("Make sure Start Time, table, and flight level coordinates are set properly")
                end

                namelistText(end + 1:end+5) = [...
                    " air_speed = 120. ! meters per second", ...
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

                namelistText(end + 1) = sprintf(" scanning_table                = './output/scanlists/%s.txt'",lower(node.Parent.Text(1:3)));
                namelistText(end + 1) = " pulse_repetition_frequency    = "+app.pulse_repetition_frequency.Value;
                namelistText(end + 1) = " pulses_per_pulse_set          = "+app.pulses_per_pulse_set.Value;
                namelistText(end + 1) = " revisits_per_acquisition_time = "+app.revisits_per_acquisition_time.Value;
                namelistText(end + 1) = " beams_per_acquisition_time    = "+app.beams_per_acquisition_time.Value;
                namelistText(end + 1) = " skip_seconds_between_scans    = 0.0";
                namelistText(end + 1) = " meters_between_gates          = "+app.meters_between_gatesEditField.Value;
                namelistText(end + 1) = " meters_to_center_of_first_gate= "+app.meters_to_center_of_first_gate.Value;
                namelistText(end + 1) = " max_range_in_meters           = "+app.max_range_in_meters.Value;
                namelistText(end + 1) = "/";

                %% wrap up
                namelistText(end + 1) = "&config_output";
                namelistText(end + 1) = "/";

                %% write to file
                writelines(namelistText, "./output/namelists/"+lower(node.Parent.Text(1:3)+"_namelist.txt"))

            end
            
        end

        % Button pushed function: InitAOSPREButton
        function InitAOSPREButtonPushed(app, event)
            
            %% Check to make sure everything has been filled in. 
            app.Messages.Text = {''};
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
                    app.Messages.Text{errorCounter} = ['Variable "', parentsString, varsToCheck{iVar}, '" is missing.'];
                    app.Messages.BackgroundColor = [0.9, 0.8, 0.8];
                    
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
                    app.Messages.Text{end+1} = ['Variable "', parentsString, app.(pathsToCheck{iPath}).Title, '" is missing.'];
                    app.Messages.BackgroundColor = [0.9, 0.8, 0.8];
                    errorCounter = errorCounter + 1;
                end
            end

            

            if errorCounter > 0
                return
            end

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
            
            app.Messages.Text = {...
                'AOSPRE-GUI completed without errors.', ...
                'Run `bash ./output/run-aospre.sh` in terminal.', ...
                'Text has been added to system clipboard.' ...
            };
            app.Messages.BackgroundColor = 'None';
            clipboard('copy','bash ./output/run-aospre.sh')


        end

        % Close request function: UIFigure
        function UIFigureCloseRequest(app, event)
            delete("./.meta/aospregui.log", 'file')
            writelines(...
                {...
                    app.SimulationGlobPattern.Value, ...
                    app.TimeGlobPattern.Value, ...
                    num2str(app.SimulationResolutionEditField.Value), ...
                    app.AOSPREPath.Value ...
                }, ...
                './.meta/aospregui.log'...
            );

            writematrix(app.UITable.Data, './.meta/table.txt')
            delete(app)
        end

        % Button pushed function: InteractiveTargetPlannerButton
        function InteractiveTargetPlannerButtonPushed(app, event)
            %% call setTargetInteractive class
            TI = setTargetInteractive(app.environment);
            uiwait(TI.fig)
            
            %% save target(s) information
            app.Targets.map = TI.Targets.map;
            app.Targets.weight = TI.Targets.weight;
            app.Targets.dwellTime = TI.Targets.dwellTime;
            


            
    
            %% update figure once targets havebeen retrieved
            app.UpdateFigure()
            
            

            

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
            w = ncread(...
                [...            % load variable at appropriate time
                    app.FileList(app.TimeIndexEditField.Value).folder, '/', ...
                    app.FileList(app.TimeIndexEditField.Value).name ...
                ], ...
                'W' ...   %variable to load
            );
            
            wSlice = w(:,:,app.SimulationLevelEditField.Value)';
            
            FP = initializeEnvironment(FP, wSlice, app.SimulationResolutionEditField.Value);
            if strcmp(app.StartOrContinueFlight.Value, 'Continue')
                % pass starting location to FlightPlanner
                x0 = app.UITable.Data(:,[1,3]);
                x = [x0(1,1), x0(:,2)'];
                y0 = app.UITable.Data(:,[2,4]);
                y = [y0(1,1), y0(:,2)'];
                
                FP.priorPath = [x',y'];
            else
                bool = 1;
                if ~isempty(app.UITable.Data)
                    uic = uiconfirm(app.UIFigure, 'Are you sure? This will delete your current flight path.', 'Delete flight path')
                    bool = strcmp(uic, 'OK');
                end

                if bool == 1
                    app.UITable.Data = [];
                else
                    return
                end
            end


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

            FP.drawGui();

            
            uiwait(FP.fig)

            app.flightPath = FP.flightPath.roi;
            app.flightPath.Parent = app.UIAxes;
            app.UpdateTableLogic()
            app.UpdateFigure();


        end

        % Button pushed function: EditpathButton
        function EditpathButtonPushed(app, event)
            
            if ~strcmp(app.EditpathButton.Text, 'Done editing')
                x0 = app.UITable.Data(:,[1,3]);
                x = [x0(1,1), x0(:,2)'];
                y0 = app.UITable.Data(:,[2,4]);
                y = [y0(1,1), y0(:,2)'];
                
                app.flightPath = images.roi.Polyline;
                app.flightPath.Position = [x',y'];
                app.flightPath.Parent = app.UIAxes;
                app.EditpathButton.Text = 'Done editing';
                app.EditpathButton.BackgroundColor = [0.5, 1, 0.5];
            else
                app.EditpathButton.Text = 'Edit path';
                app.EditpathButton.BackgroundColor = [0.95, 0.95, 0.95];
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
            app.SimulationGlobPattern.Value = [d,'/*.nc'];
            
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 820 516];
            app.UIFigure.Name = 'MATLAB App';
            app.UIFigure.CloseRequestFcn = createCallbackFcn(app, @UIFigureCloseRequest, true);

            % Create UIAxes
            app.UIAxes = uiaxes(app.UIFigure);
            title(app.UIAxes, 'Title')
            xlabel(app.UIAxes, 'X')
            ylabel(app.UIAxes, 'Y')
            zlabel(app.UIAxes, 'Z')
            app.UIAxes.Position = [500 200 320 317];

            % Create PlotDetailsPanel
            app.PlotDetailsPanel = uipanel(app.UIFigure);
            app.PlotDetailsPanel.Title = 'Plot Details';
            app.PlotDetailsPanel.Position = [499 0 321 198];

            % Create GridLayout3
            app.GridLayout3 = uigridlayout(app.PlotDetailsPanel);
            app.GridLayout3.ColumnWidth = {120, 80, 80};
            app.GridLayout3.RowHeight = {22, 22, 22, 22, 34};
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
            app.TimeIndexLabel.Visible = 'off';
            app.TimeIndexLabel.Layout.Row = 4;
            app.TimeIndexLabel.Layout.Column = 1;
            app.TimeIndexLabel.Text = 'Time Index';

            % Create TimeIndexSlider
            app.TimeIndexSlider = uislider(app.GridLayout3);
            app.TimeIndexSlider.ValueChangedFcn = createCallbackFcn(app, @TimeIndexSliderValueChanged, true);
            app.TimeIndexSlider.Visible = 'off';
            app.TimeIndexSlider.Layout.Row = 5;
            app.TimeIndexSlider.Layout.Column = [1 2];

            % Create SimulationLevelSlider
            app.SimulationLevelSlider = uislider(app.GridLayout3);
            app.SimulationLevelSlider.Orientation = 'vertical';
            app.SimulationLevelSlider.ValueChangedFcn = createCallbackFcn(app, @SimulationLevelSliderValueChanged, true);
            app.SimulationLevelSlider.Visible = 'off';
            app.SimulationLevelSlider.Layout.Row = [1 5];
            app.SimulationLevelSlider.Layout.Column = 3;

            % Create SimulationLevelEditField
            app.SimulationLevelEditField = uieditfield(app.GridLayout3, 'numeric');
            app.SimulationLevelEditField.ValueChangedFcn = createCallbackFcn(app, @SimulationLevelEditFieldValueChanged, true);
            app.SimulationLevelEditField.Visible = 'off';
            app.SimulationLevelEditField.Layout.Row = 3;
            app.SimulationLevelEditField.Layout.Column = 2;
            app.SimulationLevelEditField.Value = 1;

            % Create SimLevelLabel
            app.SimLevelLabel = uilabel(app.GridLayout3);
            app.SimLevelLabel.HorizontalAlignment = 'right';
            app.SimLevelLabel.Visible = 'off';
            app.SimLevelLabel.Layout.Row = 3;
            app.SimLevelLabel.Layout.Column = 1;
            app.SimLevelLabel.Text = 'Sim. Level';

            % Create PlottedVariableDropDownLabel
            app.PlottedVariableDropDownLabel = uilabel(app.GridLayout3);
            app.PlottedVariableDropDownLabel.HorizontalAlignment = 'right';
            app.PlottedVariableDropDownLabel.Visible = 'off';
            app.PlottedVariableDropDownLabel.Layout.Row = 1;
            app.PlottedVariableDropDownLabel.Layout.Column = 1;
            app.PlottedVariableDropDownLabel.Text = 'Plotted Variable';

            % Create PlottedVariableDropDown
            app.PlottedVariableDropDown = uidropdown(app.GridLayout3);
            app.PlottedVariableDropDown.Items = {'Set Search Pattern'};
            app.PlottedVariableDropDown.ValueChangedFcn = createCallbackFcn(app, @PlottedVariableDropDownValueChanged, true);
            app.PlottedVariableDropDown.Visible = 'off';
            app.PlottedVariableDropDown.Layout.Row = 1;
            app.PlottedVariableDropDown.Layout.Column = 2;
            app.PlottedVariableDropDown.Value = 'Set Search Pattern';

            % Create TabGroup
            app.TabGroup = uitabgroup(app.UIFigure);
            app.TabGroup.Position = [1 3 499 514];

            % Create SimulationDetailsTab
            app.SimulationDetailsTab = uitab(app.TabGroup);
            app.SimulationDetailsTab.Title = 'Simulation Details';

            % Create FileSearchPatternEditFieldLabel
            app.FileSearchPatternEditFieldLabel = uilabel(app.SimulationDetailsTab);
            app.FileSearchPatternEditFieldLabel.HorizontalAlignment = 'right';
            app.FileSearchPatternEditFieldLabel.Position = [3 454 107 22];
            app.FileSearchPatternEditFieldLabel.Text = 'File Search Pattern';

            % Create SimulationGlobPattern
            app.SimulationGlobPattern = uieditfield(app.SimulationDetailsTab, 'text');
            app.SimulationGlobPattern.ValueChangedFcn = createCallbackFcn(app, @SimulationGlobPatternValueChanged, true);
            app.SimulationGlobPattern.Tooltip = {''};
            app.SimulationGlobPattern.Placeholder = '/path/to/simulation/directory/*.nc';
            app.SimulationGlobPattern.Position = [123 454 250 22];

            % Create FileSearchFeedback
            app.FileSearchFeedback = uilabel(app.SimulationDetailsTab);
            app.FileSearchFeedback.Position = [10 319 319 33];
            app.FileSearchFeedback.Text = '';

            % Create FileTimePatternEditFieldLabel
            app.FileTimePatternEditFieldLabel = uilabel(app.SimulationDetailsTab);
            app.FileTimePatternEditFieldLabel.HorizontalAlignment = 'right';
            app.FileTimePatternEditFieldLabel.Position = [14 355 96 22];
            app.FileTimePatternEditFieldLabel.Text = 'File Time Pattern';

            % Create TimeGlobPattern
            app.TimeGlobPattern = uieditfield(app.SimulationDetailsTab, 'text');
            app.TimeGlobPattern.ValueChangedFcn = createCallbackFcn(app, @TimeGlobPatternValueChanged, true);
            app.TimeGlobPattern.Tooltip = {'Enter search pattern for output files. I''ve done my best to guess, correct if necessary. Use the following format: yyMMddHHmmss.'; 'e.g.: ''????????dHHmmss????'''; 'or'; '''?????????sssss????*'''};
            app.TimeGlobPattern.Position = [125 354 248 23];

            % Create FileBrowserButton
            app.FileBrowserButton = uibutton(app.SimulationDetailsTab, 'push');
            app.FileBrowserButton.ButtonPushedFcn = createCallbackFcn(app, @FileBrowserButtonPushed, true);
            app.FileBrowserButton.Position = [383 453 104 23];
            app.FileBrowserButton.Text = 'File Browser';

            % Create SimulationResolutionEditField
            app.SimulationResolutionEditField = uieditfield(app.SimulationDetailsTab, 'numeric');
            app.SimulationResolutionEditField.Limits = [0 Inf];
            app.SimulationResolutionEditField.ValueChangedFcn = createCallbackFcn(app, @SimulationResolutionEditFieldValueChanged, true);
            app.SimulationResolutionEditField.Placeholder = '500';
            app.SimulationResolutionEditField.Position = [125 404 248 22];
            app.SimulationResolutionEditField.Value = 1;

            % Create SimResolutionLabel
            app.SimResolutionLabel = uilabel(app.SimulationDetailsTab);
            app.SimResolutionLabel.HorizontalAlignment = 'right';
            app.SimResolutionLabel.Position = [-10 404 120 22];
            app.SimResolutionLabel.Text = 'Sim. Resolution (m)';

            % Create FlightPlanningTab
            app.FlightPlanningTab = uitab(app.TabGroup);
            app.FlightPlanningTab.Title = 'Flight Planning';

            % Create InteractiveTargetPlannerButton
            app.InteractiveTargetPlannerButton = uibutton(app.FlightPlanningTab, 'push');
            app.InteractiveTargetPlannerButton.ButtonPushedFcn = createCallbackFcn(app, @InteractiveTargetPlannerButtonPushed, true);
            app.InteractiveTargetPlannerButton.Position = [14 439 148 23];
            app.InteractiveTargetPlannerButton.Text = 'Interactive Target Planner';

            % Create ManuallydrawpathButton
            app.ManuallydrawpathButton = uibutton(app.FlightPlanningTab, 'push');
            app.ManuallydrawpathButton.ButtonPushedFcn = createCallbackFcn(app, @ExtractWaypoints, true);
            app.ManuallydrawpathButton.Tooltip = {'Select waypoints in the plot and use plot route to save them to the table. this will also plot the route in the figure.'};
            app.ManuallydrawpathButton.Position = [7 314 215 23];
            app.ManuallydrawpathButton.Text = 'Manually draw path';

            % Create UITable
            app.UITable = uitable(app.FlightPlanningTab);
            app.UITable.BackgroundColor = [1 1 1;1 1 1;0.902 0.902 0.902;0.902 0.902 0.902;1 1 1;1 1 1];
            app.UITable.ColumnName = {'Start X'; 'Start Y'; 'End X'; 'End Y'; 'Leg Time'; 'Heading'};
            app.UITable.ColumnWidth = {'fit'};
            app.UITable.RowName = {};
            app.UITable.ColumnEditable = true;
            app.UITable.CellEditCallback = createCallbackFcn(app, @UITableCellEdit, true);
            app.UITable.Tooltip = {'right click to open context menu: clear, refresh, or add new leg.'};
            app.UITable.Position = [1 194 498 107];

            % Create InteractiveFlightPlannerButton
            app.InteractiveFlightPlannerButton = uibutton(app.FlightPlanningTab, 'push');
            app.InteractiveFlightPlannerButton.ButtonPushedFcn = createCallbackFcn(app, @InteractiveFlightPlannerButtonPushed, true);
            app.InteractiveFlightPlannerButton.Position = [219 376 148 23];
            app.InteractiveFlightPlannerButton.Text = 'Interactive Flight Planner';

            % Create StartOrContinueFlight
            app.StartOrContinueFlight = uidropdown(app.FlightPlanningTab);
            app.StartOrContinueFlight.Items = {'Start new flight', 'Continue'};
            app.StartOrContinueFlight.Tooltip = {'Start planning a new flight (from scratch) or build from the end of the current path.'};
            app.StartOrContinueFlight.Position = [14 376 148 22];
            app.StartOrContinueFlight.Value = 'Start new flight';

            % Create EditpathButton
            app.EditpathButton = uibutton(app.FlightPlanningTab, 'push');
            app.EditpathButton.ButtonPushedFcn = createCallbackFcn(app, @EditpathButtonPushed, true);
            app.EditpathButton.Tooltip = {'Select waypoints in the plot and use plot route to save them to the table. this will also plot the route in the figure.'};
            app.EditpathButton.Position = [245 314 215 23];
            app.EditpathButton.Text = 'Edit path';

            % Create AOSPREDetailsTab
            app.AOSPREDetailsTab = uitab(app.TabGroup);
            app.AOSPREDetailsTab.Title = 'AOSPRE Details';

            % Create GridLayout4
            app.GridLayout4 = uigridlayout(app.AOSPREDetailsTab);
            app.GridLayout4.ColumnWidth = {59, 60, 55, 104, 153, '1x'};
            app.GridLayout4.RowHeight = {30, 30, 50, 0, 60, 150, 30, 30};
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
            app.StartTimeEditField.Layout.Column = [3 5];

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

            % Create InitAOSPREButton
            app.InitAOSPREButton = uibutton(app.GridLayout4, 'push');
            app.InitAOSPREButton.ButtonPushedFcn = createCallbackFcn(app, @InitAOSPREButtonPushed, true);
            app.InitAOSPREButton.Layout.Row = 8;
            app.InitAOSPREButton.Layout.Column = [3 5];
            app.InitAOSPREButton.Text = 'Init. AOSPRE';

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

            % Create AdvancedOptionsTab
            app.AdvancedOptionsTab = uitab(app.TabGroup);
            app.AdvancedOptionsTab.Title = 'Advanced Options';

            % Create outputformatstringLabel
            app.outputformatstringLabel = uilabel(app.AdvancedOptionsTab);
            app.outputformatstringLabel.HorizontalAlignment = 'right';
            app.outputformatstringLabel.WordWrap = 'on';
            app.outputformatstringLabel.Position = [4 453 176 32];
            app.outputformatstringLabel.Text = 'output_filename_format_string = ';

            % Create output_filename_format_string
            app.output_filename_format_string = uieditfield(app.AdvancedOptionsTab, 'text');
            app.output_filename_format_string.Tooltip = {'The naming pattern to use on the outputted CFradial netCDF files. "A" is the placeholder to AOSPRE to input time details.'; ''; 'PNL (SCN) will be replaced with the panel (scan) identifier.'};
            app.output_filename_format_string.Position = [149 434 246 22];
            app.output_filename_format_string.Value = '''("./output/PNL_SCN_",A,"_to_",A,".nc")''';

            % Create flight_level_coordinateEditFieldLabel
            app.flight_level_coordinateEditFieldLabel = uilabel(app.AdvancedOptionsTab);
            app.flight_level_coordinateEditFieldLabel.HorizontalAlignment = 'right';
            app.flight_level_coordinateEditFieldLabel.WordWrap = 'on';
            app.flight_level_coordinateEditFieldLabel.Position = [5 395 176 32];
            app.flight_level_coordinateEditFieldLabel.Text = 'flight_level_coordinate = ';

            % Create flight_level_coordinate
            app.flight_level_coordinate = uieditfield(app.AdvancedOptionsTab, 'text');
            app.flight_level_coordinate.Tooltip = {'The naming pattern to use on the outputted CFradial netCDF files. "A" is the placeholder to AOSPRE to input time details.'; ''; 'PNL (SCN) will be replaced with the panel (scan) identifier.'};
            app.flight_level_coordinate.Position = [150 376 246 22];
            app.flight_level_coordinate.Value = '"Z"';

            % Create flight_level_waypoints_vertEditFieldLabel
            app.flight_level_waypoints_vertEditFieldLabel = uilabel(app.AdvancedOptionsTab);
            app.flight_level_waypoints_vertEditFieldLabel.HorizontalAlignment = 'right';
            app.flight_level_waypoints_vertEditFieldLabel.WordWrap = 'on';
            app.flight_level_waypoints_vertEditFieldLabel.Position = [4 338 176 32];
            app.flight_level_waypoints_vertEditFieldLabel.Text = 'flight_level_waypoints_vert = ';

            % Create flight_level_waypoints_vert
            app.flight_level_waypoints_vert = uieditfield(app.AdvancedOptionsTab, 'text');
            app.flight_level_waypoints_vert.Tooltip = {'The naming pattern to use on the outputted CFradial netCDF files. "A" is the placeholder to AOSPRE to input time details.'; ''; 'PNL (SCN) will be replaced with the panel (scan) identifier.'};
            app.flight_level_waypoints_vert.Position = [149 319 246 22];
            app.flight_level_waypoints_vert.Value = '2000';

            % Create ScanningOptionsPanel
            app.ScanningOptionsPanel = uipanel(app.AdvancedOptionsTab);
            app.ScanningOptionsPanel.Tooltip = {'Some advanced '};
            app.ScanningOptionsPanel.Title = 'Scanning Options';
            app.ScanningOptionsPanel.Position = [19 33 364 266];

            % Create ScanningOptionsGrid
            app.ScanningOptionsGrid = uigridlayout(app.ScanningOptionsPanel);
            app.ScanningOptionsGrid.ColumnWidth = {39, 98, 36, '1x'};
            app.ScanningOptionsGrid.RowHeight = {23, 22, 22, 22, 22, 22, 22, 22};
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

            % Create meters_between_gatesEditField
            app.meters_between_gatesEditField = uieditfield(app.ScanningOptionsGrid, 'text');
            app.meters_between_gatesEditField.HorizontalAlignment = 'right';
            app.meters_between_gatesEditField.Layout.Row = 6;
            app.meters_between_gatesEditField.Layout.Column = 4;
            app.meters_between_gatesEditField.Value = '75';

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
            runStartupFcn(app, @startupFcn)

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