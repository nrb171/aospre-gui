function [ax, varargout] = plotEarth(ax, varargin)
    p = inputParser;    
    p.addParameter('Topo', false, @islogical);          % if true, will plot topography with green-yellow-grown colormap
    p.addParameter('Rasterize', false, @islogical);     % if you are using this same axes for many images, rasterizing will save significant time by pre-baking the contours.
    p.parse(varargin{:});
    warning off
    hold(ax,'on')

    % landareas = shaperead('landareas.shp');
    rivers = shaperead('worldrivers');
    lakes = shaperead('worldlakes');

    fid = fopen('./data/ETOPO5.DOS') ;   
    topo = fread(fid, 'bit16');
    fclose(fid);
    % X = 
    topo = reshape(topo, 360*12, []);
    topo2 = topo;
    topo2(1:180*12,:) = topo(180*12+1:end,:);
    topo2(180*12+1:end,:) = topo(1:180*12,:);
    topo = topo2;
    clear topo2
    
    %get coordinates
    % laX = {landareas.X};
    % laY = {landareas.Y};
    
    lakesX = {lakes.X};
    lakesY = {lakes.Y};
    
    riX = {rivers.X};
    riY = {rivers.Y};
    %plot them
    X = -180:0.5:179.5;
    Y = -90:0.5:89.5;
    [YY,XX] = meshgrid(Y, X);
    
    % ps = reduceShape(laX, laY);
    % plot(ax,ps, 'FaceColor','#75A44E', 'FaceAlpha',0.5)
    scale = 0.3;
    X1 = linspace(-180, 180, 360*12*scale);
    Y1 = linspace(-90, 90, 180*12*scale);
    

    M = contourc(X1,Y1,imresize(flipud(topo'), scale), [0,0]);
    plotM(ax,M,[0,0.3725,0.2706], 'LineWidth', 0.5, 'FaceColor', 'none')

    ps = reduceShape(lakesX, lakesY);
    plot(ax,ps, 'FaceColor','#6b93d6', 'FaceAlpha',0.5)

    % ps = reduceShape(riX, riY);
    % plot(ax,ps)


    % contourf(ax,ps, 'FaceColor','#61714D', 'FaceAlpha',0.5)
    % for i = 1:size(lakesX,2)
    %     ps=polyshape(lakesX{i},lakesY{i});
    %     
    % end
    for i = 1:size(riX,2)
        ps=polyshape(riX{i},riY{i});
        plot(ax,riX{i},riY{i}, 'Color','#4f4cb0')
    end
    warning on

 


    if p.Results.Topo
        plotM(ax,M,[0,0.3725,0.2706], 'LineWidth', 0.5)
        % b(b < 0) = NaN;
        % topomap1 = load('topography.mat', 'topomap1')
        % topomap1 = topomap1.topomap1;
        
        cm = demcmap([-0, 3200], 32);
        cm = [cm;flipud(cbrewer('seq', 'YlOrBr', 32))];
        % pcolor(flipud(imresize(b, 0.2)'));

        for lvl = 300:300:6000
            M = contourc(X1,Y1,imresize(flipud(topo'), scale), [lvl, lvl]);
            col = interp1(linspace(0, 6000, size(cm,1)), cm, lvl);
            
            plotM(ax,M,col, 'LineStyle', 'None')
            
        end

        % [cn,a] = contour(ax, X1, Y1,flipud(imresize(topo, 0.3)'), -0:250:3000);
        % clim([-0, 3000])
        % shading flat
        % print2(gcf, '.temp.png', 'Quality', '-r300')
    end

    

    if p.Results.Rasterize
        % rasterize the image and save to M x N x 3 matrix, front loads the workload
        fi = get(ax, 'Parent');
        oldFigUnits = fi.Units;
        oldFigPosition = fi.Position;
        olfFigColor = fi.Color;
        oldAxUnits = ax.Units;
        oldAxPosition = ax.Position;
        xlim(ax, [-180, 180])
        ylim(ax, [-90, 90])
        
        % 10 pixels per degree (~11 km)
        fi.Units = 'pixels';
        fi.Position = [0,0,360*6,180*6];
        ax.Units = 'pixels';
        ax.Position = fi.Position;
        fi.Color = [1,1,1];
        
        axis(ax, 'off')

        frame = getframe(ax);
        
        % the background is off-white, so this fixes it.
        frame.cdata(frame.cdata(:,:,1) == 240 & frame.cdata(:,:,2) == 240 & frame.cdata(:,:,3) == 240) = 255;
        hold(ax, 'off')
        h = image(ax,[-180, 180], [90, -90], frame.cdata)
        % h = rasterizedImage(h);
        varargout{1} = h;
        ax.YDir = 'normal';
        
        fi.Units = oldFigUnits;
        fi.Position = oldFigPosition;
        fi.Color = olfFigColor;
        ax.Units = oldAxUnits;
        ax.Position = oldAxPosition;
        hold(ax, 'on')
        axis(ax, 'on')
    end

    %
    xticks(ax, -180:10:180)
    xticklabels(ax, cellfun(@(x) num2str(x)+"^\circ", num2cell(-180:10:180), 'UniformOutput', false))
    yticks(ax, -90:10:90)
    yticklabels(ax, cellfun(@(x) num2str(x)+"^\circ", num2cell(-90:10:90), 'UniformOutput', false))
    grid(ax, 'on')




    function plotM(ax,M,col, varargin)
        i = 1;
        while true
            N = M(2,i); % Number of points in this line
            X = M(1,i+1:i+N); % Extract the line
            Y = M(2,i+1:i+N);
            X(end+1) = X(1); % Close the polygon
            Y(end+1) = Y(1);
            
            fill(ax, X, Y, col, varargin{:})
            i = i + N + 1; % Skip to the next line
            if i >= size(M,2)
                break
            end
        end
    end
    function ps = reduceShape(shapeX, shapeY)
        %% simplify the polygons
        for i = 1:size(shapeX,2)
            try
                X = shapeX{i}';
                Y = shapeY{i}';
                nanMask = ~isnan(X) & ~isnan(Y);
                ps2 = polyshape(X,Y);
                
                ps2Holes = ps2.holes;
                ps2Boundary = ps2.rmholes;
    
                % simplify boundary
                pout = regions(ps2Boundary);
                for j = 1:numel(pout)
                    X = pout(j).simplify.Vertices(:,1);
                    Y = pout(j).simplify.Vertices(:,2);
                    nanMask = ~isnan(X) & ~isnan(Y);
                    pReduced = reducepoly([X(nanMask), Y(nanMask)], 0.001);
                    if j == 1
                        psTemp = polyshape(pReduced(:,1), pReduced(:,2));
                    else
                        psTemp = union(psTemp, polyshape(pReduced(:,1), pReduced(:,2)));
                    end
                end
                
                
                
                % simplify holes
                for j = 1:numel(ps2Holes)
                    hReduced = reducepoly(ps2Holes(j).simplify.Vertices, 0.001);
                    % ps2Holes(j) = polyshape(hReduced(:,1), hReduced(:,2));
                    psTemp = addboundary(psTemp, hReduced(:,1), hReduced(:,2));
                end
    
                if i == 1
                    ps = psTemp;
                else
                    ps = union(ps, psTemp);
                end
            catch
                keyboard
            end
            
        end
    end
end

