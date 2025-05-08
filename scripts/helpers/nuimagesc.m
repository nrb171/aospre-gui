function [imsc] = nuimagesc(ax,xlat, xlon, var)
%NUIMAGESC - non-uniform imagesc plotting.
%   handles image plotting, but with non-orthogonal grids, such as often is
%   found in numerical weather simulations. This is substantially faster
%   than alternatives like PCOLOR(), which uses surfaces instead of images.
    xlat2 = xlat;
    xlon2 = xlon;
    
    var2 = var;
    indsi = round(linspace(1,size(xlon,1),31));
    indsj = round(linspace(1,size(xlon,2),31));
    for ii = 1:numel(indsi)-1
        for jj = 1:numel(indsj)-1
            if ii == 1 & jj == 1
                hold(ax, 'off')
            else
                hold(ax, 'on')
            end
            xlat = xlat2(...
                max(indsi(ii)-2,1):min(indsi(ii+1)+2, size(xlat2,1)), ...
                max(indsj(jj)-2,1):min(indsj(jj+1)+2, size(xlat2,2)));
            xlon = xlon2(...
                max(indsi(ii)-2,1):min(indsi(ii+1)+2, size(xlat2,1)), ...
                max(indsj(jj)-2,1):min(indsj(jj+1)+2, size(xlat2,2)));
            var = var2(...
                max(indsi(ii)-2,1):min(indsi(ii+1)+2, size(xlat2,1)), ...
                max(indsj(jj)-2,1):min(indsj(jj+1)+2, size(xlat2,2)));
    
            [py,px] = meshgrid(linspace(0,1,size(xlat,2)), linspace(0,1,size(xlat,1)));
            xRescaled = normalize(xlon(:), 1, 'range', [0,1]);
            xRescaled = reshape(xRescaled, size(xlon));
            
            
            yRescaled = normalize(xlat(:), 1, 'range', [0, 1]);
            yRescaled = reshape(yRescaled, size(xlon));
            
            d = cat(3, (yRescaled-py)*size(xRescaled,2), (xRescaled-px)*size(xRescaled,1));
             
            im = imwarp(var,-d*1, 'cubic');
            %im = im(sum(im,2) ~= 0, sum(im,1) ~= 0);
            
            alphaData = ones(size(im))*1;
            alphaData(im == 0) = 0;
            alphaData(isnan(im)) = 0;

            im(im == 0) = NaN;
            
            
            % close all
            imagesc(ax,...
                [min(xlon(:)), max(xlon(:))], ...
                [min(xlat(:)), max(xlat(:))], ...
                im', ...
                "AlphaData", alphaData' ...
            );
            
        end
    end
    
    set(ax, 'YDir','normal') 
end

