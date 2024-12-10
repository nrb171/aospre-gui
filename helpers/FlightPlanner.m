classdef FlightPlanner < handle
    properties
        Targets = struct(...
            'map',[], ...
            'weight',[], ...
            'dwellTime', []);   % Target structure spatiotemporal location
        Environment     % Environment structure
        flights         % graph of possible flights
        timeStep = 120   % time step for flight simulation (seconds)
        searchLength = 120*120    % length of search (meters)
        Aircraft = struct(...               % Aircraft structure (speed, maxRoll, maxTurbulence, Radars)
            'speed', 120, ...               % m/s
            'maxRoll', pi/3, ...            % degrees
            'maxTurbulence', 3, ...         % (m/s) range of vertical velocity over time
            'Radars', struct(...            % radar objects
                'lhs', 90, 'rhs', -90), ...   
            'optimalDistance', 5000, ...    % meters from target
            'distanceForgiveness', 4000, ...% meters from target
            'maxDistance', 45000)        % meters from target
        Costs = struct(...
            'headingOffset', 3, ...     % cost of aircraft turning
            'turbulence', 5, ...        % cost of turbulence
            'base', 16, ...             % cost of heading offset and turbulence
            'lengthDiscount', 5)        % discount for long legs         
        TrainingOptions = struct(...    % training options
            'learningRate', 0.008, ...      % rate of learning for valid mask
            'nLoops', 100, ...              % number of loops
            'nParallel', 100, ...           % number of sibling paths to calculate before updating the valid mask
            'trimStart', 0.3, ...           % fraction of iterations of nParallel to wait before trimming
            'trimRate', 4, ...             % rate to trim lowest percentile of trees after trimStart
            'trimAggression', 50, ...       % highest percentile quality of trees to keep
            'trimLearningRate', 0, ...    % rate that percentile of triming goes up (per loop above trimStart)
            'trimLength', 1, ...           % minimum length of bad trees to cut
            'desiredPathLength', 1, ...     % desired path length (hours)
            'condensationRate', 10, ...      % rate to condense network
            'condensationRadius', 4)         % radius to perform search (pixels)
        validMask       % weighted mask of valid locations in the environment


    end

    %% properties for UIFigure
    properties (Access = public)
        fig
        ax

        % training/init options
        startPointBtn
        finalizeStartBtn
        trainBtn
        stopTrainingBtn
        stopTraining = 0;

        learningRateEF
        nLoopsEF
        nParallelEF
        desiredPathLengthEF

        % progress bar
        progressBarAxes
        progressBarText

        % start line/point
        startLine



    end

    properties (Access = private)
        g = 9.81            % gravity constant
        searchExtender = 5  % multiplier for search length
        loopNumber = 0      % number of loops
        nodesInPath
        

        
        
    end
    
    %% Private methods for model training.
    methods (Access = private)
         function [costs,heading,targetsObserved] = costOfLegs(obj, oldWaypointIDs, newWaypointIDs)
                % loop through each possible node and calculate the cost of the leg
                %obj.Targets.dwellTime, observedTime, weight
                % distance from target
                % angular separation
                % turbulence
                % discount for length of leg

                % output
                % costs (matrix; oldID, newID): cost of each leg pairing
                % heading (matrix; oldID, newID): heading of each leg
                % targetsObserved(cell matrix; oldID, newID): targets observed during each leg
        
                
                costs = inf*ones([numel(oldWaypointIDs), numel(newWaypointIDs)]);
                heading = NaN(numel(oldWaypointIDs), numel(newWaypointIDs));
                targetsObserved = cell(numel(oldWaypointIDs), numel(newWaypointIDs));
                newWaypoints = obj.flights.Nodes(newWaypointIDs,:);
                oldWaypoints = obj.flights.Nodes(oldWaypointIDs,:);
                for iwpts = 1:numel(newWaypointIDs)
                    
                    try
                        newWaypoint = newWaypoints(iwpts,:);
                    catch
                        keyboard
                    end
                    for iff = 1:length(oldWaypointIDs)
                        oldWaypoint = oldWaypoints(iff,:);

                        %% calculate cost for new heading
                        distanceOfTrip = round(sqrt((newWaypoint.row - oldWaypoint.row)^2 + (newWaypoint.col - oldWaypoint.col)^2));
                        timeOfTrip = distanceOfTrip*obj.Environment.resolution/obj.Aircraft.speed;
                        heading(iff,iwpts) = atan2d(...
                            newWaypoint.row - oldWaypoint.row, ...
                            newWaypoint.col - oldWaypoint.col);
                        previousHeading = oldWaypoint.heading;
                        % r = obj.Aircraft.speed^2/(obj.g*tand(obj.Aircraft.maxRoll));
                        % theta = (obj.Aircraft.speed*obj.timeStep)/r;
                        % phi = abs(atan2d(1-cos(theta), sin(theta)));
                        headingRateOfChange = abs(wrapTo180(heading(iff,iwpts) - previousHeading))/timeOfTrip*pi/180;
                        maxHeadingRateOfChange = (obj.Aircraft.maxRoll)/obj.timeStep;
                        headingFraction = headingRateOfChange/maxHeadingRateOfChange;
                     
                        % keyboard
                        if headingFraction > 1
                            cost.heading(iff,iwpts) = inf;
                            continue
                        else
                            
                            cost.heading(iff,iwpts) = sin(headingFraction * pi/2)*obj.Costs.headingOffset;
                        end
                        
                        
                        %% calculate turbulence
                        
                        rowsOnRoute = round(linspace(oldWaypoint.row, newWaypoint.row, distanceOfTrip));
                        colsOnRoute = round(linspace(oldWaypoint.col, newWaypoint.col, distanceOfTrip));
                        indsOnRoute = sub2ind(size(obj.Environment.w), rowsOnRoute, colsOnRoute);
                        wAlongRoute = obj.Environment.w(indsOnRoute);
                        turbulenceFraction = (max(wAlongRoute) - min(wAlongRoute))/obj.Aircraft.maxTurbulence;
                        if numel((wAlongRoute)) ~= 0
                            if turbulenceFraction > 1 
                                cost.turbulence(iff,iwpts) = inf;
                                continue
                            else
                                cost.turbulence(iff,iwpts) = sin(turbulenceFraction*pi/2)*obj.Costs.turbulence;
                            end
                        else
                            cost.turbulence(iff,iwpts) = inf;
                            continue
                        end

                        %% calculate discount for observed targets.
                        radarFields = fieldnames(obj.Aircraft.Radars);
                        targetsObserved = {};
                        totalDiscount = 0;
                        targetsObserved{iff,iwpts} = [];
                        for j = 1:numel(radarFields)
                            %considering angle from radar, what targets are observed?
                            radarAngle = obj.Aircraft.Radars.(radarFields{j});
                            range = 1:1:obj.Aircraft.maxDistance/obj.Environment.resolution;
                            beam.row = round(newWaypoint.row + range*sind(radarAngle+heading(iff,iwpts))); 
                            beam.col = round(newWaypoint.col + range*cosd(radarAngle+heading(iff,iwpts)));
                            
                            validInds = beam.row > 0 & beam.row < size(obj.Environment.w,1) & beam.col > 0 & beam.col < size(obj.Environment.w,2);
                            beam.row = beam.row(validInds);
                            beam.col = beam.col(validInds);
                            indsBeam = sub2ind(size(obj.Environment.w), beam.row, beam.col);
                            beam.observed = obj.Targets.map(indsBeam);

                            
                            % pcolor(obj.Targets.map);
                            % hold on
                            % plot(beam.col, beam.row, 'r');
                            % shading flat;
                            % colorbar;
                            % print2(gcf, ".temp3.png");
                            observedTargets = unique(beam.observed(beam.observed > 0));
                            for tID = 1:nnz(observedTargets)
                                targetID = observedTargets(tID);
                                distance = find(beam.observed == targetID, 1)*obj.Environment.resolution;
                                discount = obj.Targets.weight(targetID);

                                %% calculate distance discount adjustment
                                if distance < obj.Aircraft.optimalDistance - obj.Aircraft.distanceForgiveness
                                    % if distance from 0 to min observable range
                                    distanceFraction = abs(distance - (obj.Aircraft.optimalDistance - obj.Aircraft.distanceForgiveness));
                                    distanceFraction = distanceFraction/(obj.Aircraft.optimalDistance - obj.Aircraft.distanceForgiveness);
                                    
                                elseif distance > obj.Aircraft.optimalDistance + obj.Aircraft.distanceForgiveness
                                    % if distance from max observable range to max distance
                                    distanceFraction = abs(distance - (obj.Aircraft.optimalDistance + obj.Aircraft.distanceForgiveness));
                                    distanceFraction = distanceFraction/(obj.Aircraft.maxDistance - (obj.Aircraft.optimalDistance + obj.Aircraft.distanceForgiveness));
                                else
                                    % if distance is within optimal range
                                    distanceFraction = 0;
                                end
                                discountAdjustment = sin(min(abs(distanceFraction),1)*pi/2)*discount/2;
                                
                                %% calculate dwell time discount adjustment
                                
                                newObservedTime = oldWaypoint.targetObservedTime{1}(targetID) + distanceOfTrip*obj.Environment.resolution/obj.Aircraft.speed;
                                if newObservedTime > obj.Targets.dwellTime(targetID)
                                    dwellTimeFraction = abs(newObservedTime - obj.Targets.dwellTime(targetID))/obj.Targets.dwellTime(targetID);
                                    discountAdjustment = discountAdjustment + sin(min(abs(dwellTimeFraction),1)*pi/2)*discount/2;
                                end
                                
                                totalDiscount = totalDiscount +discount - discountAdjustment;
                                try
                                    targetsObserved{iff,iwpts}(end+1) = targetID;
                                catch ME
                                    keyboard
                                end
                            end
                        end
                        cost.targetDiscount(iff,iwpts) = -totalDiscount;

                        %% discount for length
                        p = shortestpath(obj.flights, 1, oldWaypointIDs(iff));
                        if isempty(p)
                            cost.lengthDiscount(iff, iwpts) = 0;
                        else
                            cost.lengthDiscount(iff, iwpts) = -obj.Costs.lengthDiscount*sin(min(numel(p)/obj.nodesInPath, 1)*pi/2);
                        end

                        if isempty(targetsObserved{iff, iwpts})
                            cost.base(iff,iwpts) = obj.Costs.base;
                        else
                            cost.base(iff,iwpts) = 0;
                        end
                        
                        %% calculate total cost
                        costFields = fieldnames(cost);
                        costsTemp = 0;
                        for j = 1:numel(costFields)
                            costsTemp = costsTemp + cost.(costFields{j})(iff,iwpts);
                        end
                        
                        costs(iff, iwpts) = costsTemp;
                        % keyboard
                    end 
                end    
            end

        function obj = addNewEdge(obj, newParentID, newChildID, cost, targetsObserved)
            % When ading a new edge, recalculate the metadata of the child node
            % this function automatically updates the graph.
            

            % find old parent/edge from new child node
            oldEdgeID = find(obj.flights.Edges.EndNodes(:,2) == newChildID);
            newGrandparentID = obj.flights.Nodes(oldEdgeID,:).parent;
            if ~isempty(oldEdgeID)
                obj.flights = rmedge(obj.flights, oldEdgeID);
            end

            % add new edge
            obj.flights = addedge(obj.flights, newParentID, newChildID, cost);
            obj.flights.Nodes(newChildID,:).parent = newParentID;

            % add heading
            childRow = obj.flights.Nodes(newChildID,:).row;
            childCol = obj.flights.Nodes(newChildID,:).col;

            parentRow = obj.flights.Nodes(newParentID,:).row;
            parentCol = obj.flights.Nodes(newParentID,:).col;
            heading = atan2d(childRow - parentRow, childCol - parentCol);
            obj.flights.Nodes(newChildID,:).heading = heading;

            % calculate path from origin to new node
            pathToOrigin = shortestpath(obj.flights, 1, newChildID);

            

            % calculate distance from origin to new node
            try
                distanceOfNewEdge = sqrt(...
                    (obj.flights.Nodes(newChildID,:).row - obj.flights.Nodes(pathToOrigin(end-1),:).row)^2 + ...
                    (obj.flights.Nodes(newChildID,:).col - obj.flights.Nodes(pathToOrigin(end-1),:).col)^2) * obj.Environment.resolution;
            catch ME

            end
        
            obj.flights.Nodes(newChildID,:).totalDistance = ...
                obj.flights.Nodes(pathToOrigin(end-1),:).totalDistance + distanceOfNewEdge;

            % calculate total elapsed time
            obj.flights.Nodes(newChildID,:).elapsedTime = ...
                obj.flights.Nodes(pathToOrigin(end - 1),:).totalDistance/obj.Aircraft.speed + distanceOfNewEdge/obj.Aircraft.speed;
            
            % calculate target observed time
            for iff = 1:numel(targetsObserved)

                targetID = targetsObserved(iff);
                obj.flights.Nodes.targetObservedTime{newChildID}(targetID) = ...
                obj.flights.Nodes.targetObservedTime{pathToOrigin(end-1)} + distanceOfNewEdge/obj.Aircraft.speed;
            end
        end
        function obj = addPath(obj)
            %% modified RRT* algorithm to find the trajectory.
            %% 1. Randomly sample a point in the environment
            warning('off')
            [~,ps] = BellmanFord(obj.flights);
            warning('on')
            obj.nodesInPath = max(cellfun(@numel, [ps{:}]));


            availableInds = find(obj.validMask);
            nParallelMod = min(obj.TrainingOptions.nParallel, 10*height(obj.flights.Nodes));
            
            randInd = randsample(availableInds, nParallelMod, true, obj.validMask(availableInds));

            % add new table for randomly sampled nodes
            [row, col] = ind2sub(size(obj.validMask), randInd);
            tab = table();
            tab.row = row;
            tab.col = col;
            tab.heading = NaN(size(row));              % updated in 4.
            tab.totalCost = NaN(size(row));            % updated in 4.
            tab.parent = NaN(size(row));               % updated in 4.

            %update search length based on the number of loops
            searchLengthHelper = obj.searchLength;%/(obj.TrainingOptions.nLoops-obj.loopNumber)*obj.TrainingOptions.nLoops;


            % path metadata
            tab.targetObservedTime = cell(size(row)); %
            for i = 1:numel(row)
                tab.targetObservedTime{i} = zeros(size(obj.flights.Nodes(1,:).targetObservedTime));
            end
                % {zeros(size(obj.flights.Nodes(1,:).targetObservedTime))};    % updated in 4.]
            tab.elapsedTime = NaN(size(row));          % updated in 4.
            tab.totalDistance = NaN(size(row));             % updated in 4.            
            
            %% 2. Find the node within searchLength that has the lowest cost that is valid
            dists = sqrt((obj.flights.Nodes.row - repmat(tab.row, [1,height(obj.flights.Nodes)])').^2 + ...
                (obj.flights.Nodes.col - repmat(tab.col, [1,height(obj.flights.Nodes)])').^2)*obj.Environment.resolution;
            [distsMin,distsMinInd]=min(dists,[],1);

            % adjust distance if distance is too from potential parent nodes
            parentRow = obj.flights.Nodes(distsMinInd,:).row;
            parentCol = obj.flights.Nodes(distsMinInd,:).col;
            headingTemp = atan2d(row - parentRow, col - parentCol);
            newRow = round(parentRow + min(distsMin',searchLengthHelper).*sind(headingTemp)/obj.Environment.resolution);
            newCol = round(parentCol + min(distsMin',searchLengthHelper).*cosd(headingTemp)/obj.Environment.resolution);
            tab.row = newRow;
            tab.col = newCol;

            % recalculate dists
            dists = sqrt((obj.flights.Nodes.row - repmat(newRow, [1,height(obj.flights.Nodes)])').^2 + ...
                (obj.flights.Nodes.col - repmat(newCol, [1,height(obj.flights.Nodes)])').^2)*obj.Environment.resolution;
        
            %% 3. calculate cost of leg                
            obj.flights = addnode(obj.flights, tab);
            newWaypointIDs = height(obj.flights.Nodes)-nParallelMod+1:height(obj.flights.Nodes);
            nodesToBeRemoved = [];
            priorParentID = [];
            for ii=1:nParallelMod

                newWaypointID = newWaypointIDs(ii);
                oldWaypointIDs = find(dists(:,ii) < searchLengthHelper*1.2);

                % calculate costs
                [costs, ~, targetsObserved] = costOfLegs(obj, oldWaypointIDs, newWaypointID);
            
            %% 4. choose lowest cost parent node, add edge from parent to new node
                [lowestCost, lowCostInd] = min(costs(:));

                % if no valid path is found, remove the last node and return
                if isinf(lowestCost)
                    nodesToBeRemoved(end+1) = newWaypointID;
                    continue
                end

                % update the new node with metadata from the parent node
                
                [parentRow, ~] = ind2sub(size(costs), lowCostInd);
                tab.parent(ii) = oldWaypointIDs(parentRow);


        
                % add edge from parent to new node    
                if isnan(lowestCost)
                    error('reset network and try again')
                end
                obj = addNewEdge(obj, oldWaypointIDs(parentRow), newWaypointID, lowestCost, targetsObserved{lowCostInd});
            
            %% 5. If adding children to lowest cost node is valid, check other nodes within search radius to add as potential children
                clear costs targetsObserved
                newChildren = oldWaypointIDs(oldWaypointIDs ~= oldWaypointIDs(lowCostInd));
                if ~isempty(newChildren)
                    for i = 1:length(newChildren)
                        
                        
                        
                        [costs, ~, targetsObserved] = costOfLegs(obj, newWaypointID, newChildren(i));
                        
                        
                        oldEdgeID = find(obj.flights.Edges.EndNodes(:,2) == newChildren(i));
                        
                        oldEdgeCost = obj.flights.Edges.Weight(oldEdgeID);
                        if costs < oldEdgeCost
                            % check to see if the new grandchildren headings are valid with the new grandparent
                            grandChildren = successors(obj.flights, newChildren(i));
                            if ~isempty(grandChildren)
                                legacyHeadings = obj.flights.Nodes.heading(grandChildren);
                                newHeading = atan2d(obj.flights.Nodes.row(newChildren(i)) - obj.flights.Nodes.row(newWaypointID), ...
                                    obj.flights.Nodes.col(newChildren(i)) - obj.flights.Nodes.col(newWaypointID));
                                if any(abs(wrapTo180(legacyHeadings - newHeading)) > obj.Aircraft.maxRoll*180/pi)
                                    continue
                                end
                            else
                                continue
                            end

                            % record ID of parent whose child is being
                            % stolen (yes, very dramatic. I promise
                            % the nodes don't feel)
                            priorParentID(end+1) = predecessors(obj.flights, newChildren(i));
                            

                            % if valid, add new edge    
                            obj = addNewEdge(obj, newWaypointID, newChildren(i), costs, targetsObserved{1}); 
                               
                        end
                    end
                end
            end

            % if prior parent doesn't regain a child, *remove* the prior
            % parent. (wow, who knew graph theory was so violent.)
            
            for ipp = 1:numel(priorParentID)
                
                if numel(successors(obj.flights, priorParentID(ipp)))==0
                    
                    nodesToBeRemoved(end+1) = priorParentID(ipp);
                end
            end
            
            
            

            %remove nodes after edge assignment to keep indices consistent
            if ~isempty(nodesToBeRemoved)
                obj.flights = rmnode(obj.flights, nodesToBeRemoved);
            end

            

            %5.1 
            % 6. update the valid mask
            
            obj=obj.updateValidMask();
            obj.loopNumber = obj.loopNumber + 1;
            

           
        end

        function obj = updateValidMask(obj)
            
            %% update the valid mask based on the current location of the aircraft
            % inds = find(ones(size(obj.Environment.w)));
            % [rowEnv, colEnv] = ind2sub(size(obj.Environment.w), inds);

            % row = obj.flights.Nodes(end,:).row;
            % col = obj.flights.Nodes(end,:).col;
            % heading = obj.flights.Nodes(end,:).heading;

            
            % angles = wrapTo2Pi(atan2(rowEnv - row, colEnv - col) - heading*pi/180);
            % % dists = sqrt((colEnv-col).^2 + (rowEnv-row).^2)*obj.Environment.resolution;
            % angles = abs(wrapToPi(angles));
            % anglesMask = angles < obj.Aircraft.maxRoll;
            % anglesMask = reshape(anglesMask, size(obj.Environment.w));
            % % keyboard
            % obj.validMask(anglesMask & obj.validMask < 0.5) = 0.5;
            
            % fig = figure('Units','inches','Position',[0 0 3 2]);
            % contourf(obj.validMask, 'LineStyle', 'none');
            % shading flat;
            % colorbar;
            % print2(gcf, ".temp4.png");
            % close(fig);

            
            
            % calculate angles from center of arc
            %r = obj.Aircraft.speed^2/(obj.g*tand(obj.Aircraft.maxRoll));
            %theta = (obj.Aircraft.speed*obj.timeStep)/r;
            %phi = abs(atan2(1-cos(theta), sin(theta)));
            %calculate angles of arc from aircraft perspective

            
            % e = findedge(obj.flights, obj.flights.Nodes(end,:).parent, height(obj.flights.Nodes));
            % if e~=0
            %     cost = obj.flights.Edges.Weight(e);

            %     if cost < 0
                    
            %         % calculate mask of valid locations
            %         mask = zeros(size(obj.Environment.w));
            %         mask(row, col) = 1;
                    
            %         filtGauss = fspecial('disk', round(obj.searchLength/obj.Environment.resolution*3));
            %         motionFilter = fspecial('motion', obj.searchLength/obj.Environment.resolution*3, heading-90);
            %         % inds = find(~isnan(motionFilter));
            %         % [row, col] = ind2sub(size(motionFilter), inds);
            %         % angles = abs(wrapTo180(atan2d(row - size(motionFilter,1)/2, col - size(motionFilter,2)/2)-heading-180));
            %         % motionFilter(angles > 90) = 0;
            %         mask = imfilter(mask, filtGauss);
            %         mask = imfilter(mask, motionFilter);
            %         mask = single(mask>0);
            %         mask(mask~=0) = mask(mask~=0)*obj.TrainingOptions.learningRate;
            %         obj.validMask(mask~=0) = max(min(obj.validMask(mask~=0)+mask(mask~=0),1), 0);
            %         obj.validMask = imclose(obj.validMask, ones(5));
                    
            %         obj.validMask = imgaussfilt(obj.validMask, 10);
            %         %obj.validMask = imerode(obj.validMask, strel('disk', 1));
            %         % if max(obj.validMask(:)) == 1
            %         %     obj.validMask = obj.validMask - obj.TrainingOptions.learningRate/10;%! this may be turned off
                    
            %         %     obj.validMask(obj.validMask<0) = 0;
                        
            %         % end
                    
            %     end
                
            % end
            



            % new method
            % find all edges less than zero
            
            edges = find(obj.flights.Edges.Weight < min(prctile(obj.flights.Edges.Weight, 5), obj.Costs.base));
            parents = obj.flights.Edges.EndNodes(edges,1);
            children = obj.flights.Edges.EndNodes(edges,2);
            parentRows = obj.flights.Nodes.row(parents);
            parentCols = obj.flights.Nodes.col(parents);
            childRows = obj.flights.Nodes.row(children);
            childCols = obj.flights.Nodes.col(children);

            meanRows = round(mean([parentRows, childRows],2));
            meanCols = round(mean([parentCols, childCols],2));
            meanInds = sub2ind(size(obj.Environment.w), meanRows, meanCols);

            maskTemp = zeros(size(obj.Environment.w));
            maskTemp(meanInds) = 1;     
            maskTemp = double(maskTemp);
            % keyboard
            theGreatEqualizer = imfilter((maskTemp), double(fspecial('disk', round(2*obj.searchLength/obj.Environment.resolution))>0))>0;
            % maskTemp(maskTemp>0)=maskTemp(maskTemp>0)./sqrt(theGreatEqualizer(maskTemp>0));
            

        
            
            % apply a gaussian filter to the mask
            
            
            % maskTemp = imgaussfilt(maskTemp, round(obj.searchLength/obj.Environment.resolution)*2)*10;
            obj.validMask = single(ones(size(obj.Environment.w)))*0.5;
            obj.validMask(~theGreatEqualizer) = 0.5-obj.TrainingOptions.learningRate*obj.loopNumber;
            
            % obj.validMask = adapthisteq(obj.validMask);

            
            if max(obj.validMask(:)) > 1
                adjustmentFactor = max(obj.validMask(:))-1;
                obj.validMask = obj.validMask - adjustmentFactor;
            end


            obj.validMask(obj.validMask < 0) = 0.0001;
            
            

        end % updateValidMask


        function obj = condenseNetwork(obj)
            
            %% condensation subroutine
            % calculate pairwise distance from all nodes
            
            rowmat = squareform(pdist(obj.flights.Nodes.row));
            colmat = squareform(pdist(obj.flights.Nodes.col));
            distmat = (rowmat.^2 + colmat.^2).^0.5;
            [nrow, ncol] = find(distmat < obj.TrainingOptions.condensationRadius & distmat ~=0);
            
            %select nodes to optimize and loop through all of them
            nodesToWorkOn = unique(nrow);
            nodesRemoved = 0;
            nodesToRemove = [];
            alreadyWorkedOn = [];
            for ii = 1:numel(nodesToWorkOn)
                % skip if already done work
                % find all fellow parents
                parents = ncol(nrow == nodesToWorkOn(ii));
                parents(end + 1) = nodesToWorkOn(ii);

                parents(ismember(parents, alreadyWorkedOn)) = [];

                % remove from working pool
                alreadyWorkedOn = [alreadyWorkedOn(:); parents(:)];

               

                % find children of parents
                children = [];
                for jj = 1:numel(parents)
                    childrenTemp = successors(obj.flights, parents(jj));
                    children = [children; childrenTemp];

                end
                
                % calculate average lowest cost child
                [costs, ~, targetsObs] = costOfLegs(obj, parents, children);


                costsMean = mean(costs, 2);

                % skip if unable to optimize.
                if nnz(isinf(costsMean)) == numel(costsMean)
                    continue
                end

                % choose best node
                [~,favoriteParentInd] = min(costsMean);
                favoriteParent = parents(favoriteParentInd);

                % remove best node from children pool
                parents(favoriteParentInd) = [];

                % add new edges from parents to best child
                for jj = 1:numel(children)
                    obj = addNewEdge(obj, ...
                        favoriteParent, children(jj), ...
                        costs(favoriteParentInd, jj), ...
                        targetsObs{favoriteParentInd, jj});
                end

                % remove old parent nodes
                
                nodesToRemove = [nodesToRemove(:);parents(:)];
                nodesRemoved = nodesRemoved + numel(parents);

            end
            obj.flights = rmnode(obj.flights, nodesToRemove);
            fprintf("Condensed "+num2str(nodesRemoved)+" nodes.\n")

        end
        function obj = trimNetwork(obj)
            
            %% find all children of every node            
            G = obj.flights;
            GF = flipedge(G);
            d = distances(GF);
            ce = centrality(G, 'outdegree');
            
            
            warning('off')
            [dbf, pbf] = BellmanFord(G);
            warning('on')
            npbf = cellfun(@numel, [pbf{:}]);
            allCosts = dbf./npbf';

            
            
            
            cutLevel = prctile(allCosts,obj.TrainingOptions.trimAggression);
            
            % cutLevel = prctile(allCosts, obj.TrainingOptions.trimAggression);

            indsToCheck = find(ce == 0);

            nodesToBeRemoved = [];
            parentHitCount = zeros(size(ce));
            for i = 2:numel(indsToCheck)
                % calculate paths and distances along path from node to root
                p = shortestpath(GF, indsToCheck(i), 1);
                dd = gradient(d(indsToCheck(i), p));
                stopInd = find(dd < cutLevel,1)-1;

                
                
                if isempty(stopInd)
                    stopInd = numel(p)-1;
                end
                
                
                if stopInd > obj.TrainingOptions.trimLength 
                %    keyboard
                    branchPoint = find(ce(p(1:stopInd+1)) > parentHitCount(p(1:stopInd+1))+1, 1);
                    parentHitCount(p(1:min(branchPoint, numel(p)))) = parentHitCount(p(1:min(branchPoint, numel(p)))) + 1;
                    nodesToBeRemoved = [nodesToBeRemoved, p(1:branchPoint-1)];
                end

                % nodesToBeRemoved = nodesToBeRemoved(1:min(round(numel(obj.TrainingOptions.trimRate*obj.TrainingOptions.nParallel*0.85))), numel(nodesToBeRemoved));

            end

            % remove bad trees
            G = rmnode(G, nodesToBeRemoved);

            fprintf(char("Trimmed "+num2str(nnz(nodesToBeRemoved))+" nodes.\n"))

            % remove nodes that are children of bad trees (shouldn't happen
            % with new implementation)
            distFromStart = distances(G, 1);
            G = rmnode(G, find(isinf(distFromStart)));

            % save back into object
            obj.flights = G;
        end % trimGraph
    end % private methods for model training

    %% Public methods for flight planning.
    methods (Access = public)
        function obj = initializeFlight(obj, startLocation, heading)
            %% initialize the graph, the first node is the starting location.
            
            % node-specific metadata
            
            obj.flights = digraph();
            tab = table();
            tab.row = startLocation(1);
            tab.col = startLocation(2);
            tab.heading = heading;
            tab.parent = NaN;
            tab.totalCost = 10;

            % path metadata
            
            tab.targetObservedTime = {zeros([nnz(unique(obj.Targets.map)),1])};
            tab.elapsedTime = 0;
            tab.totalDistance = 0;

            % add the first node
            obj.flights = addnode(obj.flights, tab);

            
            obj.validMask = single(ones(size(obj.Environment.w))*0.5);


            %% to turn on initial weighting, use the following
            % obj.validMask = single(zeros(size(obj.Environment.w))*0.5);
            % inds = find(~obj.validMask);

            % [rowEnv, colEnv] = ind2sub(size(obj.Environment.w), inds);

            % row = obj.flights.Nodes(end,:).row;
            % col = obj.flights.Nodes(end,:).col;
            % heading = obj.flights.Nodes(end,:).heading;

            % angles = wrapTo2Pi(atan2(rowEnv - row, colEnv - col) - heading*pi/180);
            % angles = abs(wrapToPi(angles));
            % angles = angles < obj.Aircraft.maxRoll;
            % angles = reshape(angles, size(obj.Environment.w));
            % obj.validMask(angles & obj.validMask == 0.5) = 1;
            % obj.validMask(obj.validMask == 0) = 0.01;

            % fig = figure('Units','inches','Position',[0 0 3 2]);
            % contourf(obj.validMask, 'LineStyle', 'none');
            % hold on
            % plot(startLocation(2), startLocation(1), 'r*');
            % shading flat;
            % colorbar;
            % print2(gcf, ".temp4.png");
            % keyboard

        end

        function obj = initializeEnvironment(obj, verticalVelocity, resolution)
            %% 2D vertical velocity environment
            % transpose the environment to match 

            obj.Environment.resolution = resolution;
            obj.Environment.w = verticalVelocity;
            obj.Targets.map = zeros(size(obj.Environment.w));
        end
        function obj = addTarget(obj, target, varargin)
            %% pixels where the target is located
            % Weight (10): reduces the cost of the path when being observed
            % DwellTime (inf; seconds): time to spend observing the target
            %   if ObservedTime (seconds) is beyond DwellTime, reduce the weight of the target at a rate of 1-cos((ObservedTime-DwellTime)/DwellTime * pi/2)
            p = inputParser;
            
            addParameter(p, 'Weight', 10);
            addParameter(p, 'ObservedTime', 0); % time the target was observed      
            addParameter(p, 'DwellTime', inf); % time to spend at the target

            parse(p, varargin{:});
            if ~exist("obj.Targets.map", "var")
                obj.Targets.map = zeros(size(obj.Environment.w));
                targetID = 1;
            else
                targetID = unique(obj.Targets.map(:));
            end
            % keyboard
            inds = find(target);
            obj.Targets.map(inds) = targetID;
            obj.Targets.dwellTime(targetID) = p.Results.DwellTime;
            obj.Targets.observedTime(targetID) = p.Results.ObservedTime;
            obj.Targets.weight(targetID) = p.Results.Weight;
        end

        function obj = train(obj, varargin)
            % train the model to find the optimal path
            p = inputParser;
            addParameter(p, 'Verbose', false);
            addParameter(p, 'Verbosity', 10);
            
            obj = obj.refreshProgressBar(0);
            parse(p, varargin{:});


            startLoop = obj.loopNumber;
            for i = startLoop:obj.TrainingOptions.nLoops
                obj = addPath(obj);

                if any(round(obj.TrainingOptions.nLoops*obj.TrainingOptions.trimStart):obj.TrainingOptions.trimRate:obj.TrainingOptions.nLoops == i)
                    obj = trimNetwork(obj);
                end

                if mod(obj.loopNumber, obj.TrainingOptions.condensationRate) == 0
                    obj = condenseNetwork(obj);
                end

                if p.Results.Verbose == 1 && mod(i, p.Results.Verbosity) == 0
                   obj.plotNetwork();
                end
                
                

                fprintf("Training loop %d of %d\n", i, obj.TrainingOptions.nLoops);
                obj=obj.refreshProgressBar(i/obj.TrainingOptions.nLoops);
                
                if obj.stopTraining == 1
                    obj.stopTraining = 0;
                    break
                end
            end
            
           
        end

        function obj = plotNetwork(obj, varargin)
            p = inputParser;
            addParameter(p, 'SaveLocation', '.temp.png');
            parse(p, varargin{:});

            if isa(obj.fig, 'matlab.ui.Figure')
                ax = obj.ax;
                
            else
                fig = figure('Units', 'inches', 'Position', [0, 0, 4, 3]);
                ax = axes(fig);
                
            end
            hold(ax, 'off');
            contour(ax, obj.validMask, 1, 'LineColor', 'r');
            hold(ax, 'on');
            contour(ax, obj.Targets.map, 1, 'LineColor', 'magenta');
            xy1 = [obj.flights.Nodes.col, obj.flights.Nodes.row];
           
            costOfLeg = 1;
            for i = 2:height(obj.flights.Nodes)
                costOfLeg(i) = distances(obj.flights, predecessors(obj.flights, i), i);
            end
            plot(ax,obj.flights, 'XData', xy1(:,1), 'YData', xy1(:,2), 'NodeLabel', {}, 'EdgeAlpha', 0.1, 'NodeColor', 'blue', 'MarkerSize', (costOfLeg-min(costOfLeg)+0.01)/10);
            daspect(ax, [1,1,1])

            %% plot 5 shortest paths
            start = 1;
            
            totalTimeElapsed  =  obj.TrainingOptions.desiredPathLength; %seconds since start (controls length of leg)
            
            roots = find(obj.flights.Nodes.elapsedTime > totalTimeElapsed*3600);
            if isempty(roots)
                totalTimeElapsed  =  0; %seconds since start (controls length of leg)
                roots = find(obj.flights.Nodes.elapsedTime > totalTimeElapsed);
            end
            G = obj.flights;
            
            d = distances(G,1); %distances from start to all nodes.
            for i = 1:numel(roots)
                [p, ~] = shortestpath(obj.flights, start, roots(i));
                pathLengths(i) = prctile(gradient(d(p)), 50);
            end
            
            [~, idx] = sort(pathLengths, 'ascend');
            for i = 1:min(10, numel(idx))
                path = shortestpath(obj.flights, start, roots(idx(i)));
                xy = [obj.flights.Nodes.col(path), obj.flights.Nodes.row(path)];
                plot(ax, xy(:,1), xy(:,2), 'Color', [1,0.8,0.5], 'LineWidth', 2);
            end
            

            if isa(obj.fig, 'matlab.ui.Figure')
                obj.ax = ax;
            else
                print2(fig, p.Results.SaveLocation);
                close(fig);
            end
        end
    end

    %% methods for GUI
    methods (Access = private)
        function obj = drawStart(obj, option)
            if strcmp(option, 'set')
                obj.refreshFigure()
                obj.startLine = drawline(obj.ax);
                
            elseif strcmp(option, 'finalize')
                obj.initializeFlight([obj.startLine.Position(1,2), obj.startLine.Position(1,1)], atan2d(obj.startLine.Position(2,2) - obj.startLine.Position(1,2), obj.startLine.Position(2,1) - obj.startLine.Position(1,1)));
                obj.validMask = single(ones(size(obj.Environment.w))*0.5);
                obj.startLine.Parent = [];
                obj.loopNumber = 1;
                obj.refreshFigure();

            end

        end
        function obj = editFieldCallback(obj, field)
            % update the training options when the edit field is changed
            obj.TrainingOptions.(field) = obj.(field+"EF").Value;
        end
    end

    %% component creation/update
    methods (Access = private)
        function obj = createComponents(obj)
            % Draw the UI and assign logic functions to items, update
            % target map.

            %% UIFigure and UIAxes
            obj.fig = uifigure('Position',[10 10 500 275]*2);        
            obj.ax = uiaxes(obj.fig, "Position",[150, 10, 340, 255]*2);

            %% flight path buttons (clear, erase, finalize, etc.)
            % setPositionBtn = uibutton(obj.fig, ...
            %     "ButtonPushedFcn", @(src,event) obj.drawStart(), ...
            %     'Text', 'Set initial position', ...
            %     'Tooltip', "Set/reset start position", ...
            %     'Position', [150, 200, 100, 25]);
            doneBtn = uibutton(obj.fig, ...
                "ButtonPushedFcn", @(src,event) close(obj.fig), ...
                'Text', 'Done/Exit', ...
                'Tooltip', "Close the GUI.", ...
                'Position', [30, 200, 100, 25], ...
                'BackgroundColor', [1 0.5 0.5]);
            % finalizeBtn = uibutton(obj.fig, ...
            %     "ButtonPushedFcn", @(src,event) obj.finalizeTarget(), ...
            %     'Text', 'Finalize Target', ...
            %     'Tooltip',"Select desired flight path and prepare export to AOSPRE-GUI.", ...
            %     'Position', [150, 250, 100, 25]); 
                
            obj.desiredPathLengthEF = uieditfield(obj.fig, 'numeric', ...
                "Limits",[0, 1000], ...
                "Position", [150, 450, 100, 25], ...
                "Tooltip", "How long do you want the path to be (hours)? ([0,10]; "+string(obj.TrainingOptions.desiredPathLength)+"=default)", ...
                "Value", obj.TrainingOptions.desiredPathLength, ...
                "ValueChangedFcn", @(src,event) obj.editFieldCallback("desiredPathLength"));
            uilabel(obj.fig, ...
                "Text", "Desired Path Length", ...
                "Position", [30, 450, 120, 25]);
            obj.learningRateEF = uieditfield(obj.fig, 'numeric', ...
                "Limits",[0 1], ...
                "Position", [150, 400, 100, 25], ...
                "Tooltip", "How aggressively do you want to prioritize high-quality leg optimization? ([0,0.01]; "+string(obj.TrainingOptions.learningRate)+"=default)", ...
                "Value", obj.TrainingOptions.learningRate, ...
                "ValueChangedFcn", @(src,event) obj.editFieldCallback("learningRate"));
            uilabel(obj.fig, ...
                "Text", "Learning Rate", ...
                "Position", [30, 400, 120, 25]);

            obj.nLoopsEF = uieditfield(obj.fig,"numeric", ...
                "Limits",[100, 1000], ...
                "Position", [150, 350, 100, 25], ...
                "Tooltip", "How many iterations do you want to run? ([100,1000]; "+string(obj.TrainingOptions.nLoops)+"=default)", ...
                "Value", obj.TrainingOptions.nLoops, ...
                "ValueChangedFcn", @(src,event) obj.editFieldCallback("nLoops"));
            uilabel(obj.fig, ...
                "Text", "Number of Iterations", ...
                "Position", [30, 350, 120, 25]);

            obj.nParallelEF = uieditfield(obj.fig,"numeric", ...
                "Limits",[1, 1000], ...
                "Position", [150, 300, 100, 25], ...
                "Value", obj.TrainingOptions.nParallel, ...
                "Tooltip", "How many sibling paths do you want to calculate before updating the valid mask? ([1,1000]; "+string(obj.TrainingOptions.nParallel)+"=default)", ...
                "ValueChangedFcn", @(src,event) obj.editFieldCallback("nParallel"));
            uilabel(obj.fig, ...
                "Text", "Number of Parallel Paths", ...
                "Position", [30, 300, 120, 25]);
            

            %% training settings
            obj.startPointBtn = uibutton(obj.fig, ...
                "ButtonPushedFcn", @(src,event) obj.drawStart('set'), ...
                'Text', 'Set start for flight. (click and drag)', ...
                'Position', [30, 250, 200, 25]);
            
            obj.finalizeStartBtn = uibutton(obj.fig, ...
                "ButtonPushedFcn", @(src,event) obj.drawStart('finalize'), ...
                'Text', 'Finalize Start position/reset network.', ...
                'Position', [30, 200, 200, 25]);

            obj.trainBtn = uibutton(obj.fig, ...
                "ButtonPushedFcn", @(src,event) obj.trainBtnCallback(), ...
                'Text', 'Begin/continue training.', ...
                'Position', [30, 150, 200, 25]);

            obj.stopTrainingBtn = uibutton(obj.fig, ...
                "ButtonPushedFcn", @(src,event) obj.stopTrainingBtnCallback(), ...
                'Text', 'Stop training.', ...
                'Position', [30, 100, 200, 25]);
            
            
            %% Progress bar initialization
            obj.progressBarAxes = uiaxes(obj.fig, ...
                'Position', [30, 50, 200, 25], ...
                'Visible', 'off', ...
                'XLim',[0,1], ...
                'YLim',[0,1], ...
                'XTick',[], ...
                'YTick',[], ...
                'Color', 'none');

            obj.progressBarText = uilabel(obj.fig, ...
                'Text', 'Training Progress', ...
                'Position', [30, 25, 200, 25], ...
                'Visible', 'off');

        end

        function obj = refreshFigure(obj)
            hold(obj.ax, 'off');
            pcolor(obj.ax, obj.Environment.w);
            shading(obj.ax, 'flat');

            try
                hold(obj.ax, 'on');
                quiver(obj.ax, ...
                    obj.startLine.Position(1,1), obj.startLine.Position(1,2), ...
                    obj.startLine.Position(2,1) - obj.startLine.Position(1,1), ...
                    obj.startLine.Position(2,2) - obj.startLine.Position(1,2), ...
                    'Color', 'red', 'LineWidth', 2);
            catch ME
                hold(obj.ax, 'off');
            end
            
            


        end

        function obj = trainBtnCallback(obj)
            %% update GUI to reflect training status
            
            obj = obj.train('Verbose', 1,'Verbosity', 2);
        end

        function obj = stopTrainingBtnCallback(obj)
            %% update GUI to reflect training status
            obj.stopTraining = 1;
        end

        function obj = refreshProgressBar(obj, percent)
            if isa(obj.fig, 'matlab.ui.Figure')
                obj.progressBarAxes.Visible = 'on';
                obj.progressBarText.Visible = 'on';
                hold(obj.progressBarAxes, 'off');
                patch(obj.progressBarAxes, [0, percent, percent, 0], [0, 0, 1, 1], 'blue');
                obj.progressBarText.Text = string(percent*100) + "%";
                drawnow limitrate
            end
        end
    end

    %% Constructor methods
    methods (Access = public)
        function obj = drawGui(obj)
            obj.createComponents();
            obj.refreshFigure();

        end
    end
end