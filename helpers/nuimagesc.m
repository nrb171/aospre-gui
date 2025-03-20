function [imsc] = nuimagesc(ax,xlat, xlon, var)
%NUIMAGESC - non-uniform imagesc plotting.
%   handles image plotting, but with non-orthogonal grids, such as often is
%   found in numerical weather simulations. This is substantially faster
%   than alternatives like PCOLOR(), which uses surfaces instead of images.

    [py,px] = meshgrid(1:size(xlat,2), 1:size(xlat,1));
    xRescaled = normalize(xlon(:), 1, 'range', [1, size(xlon,1)]);
    xRescaled = reshape(xRescaled, size(xlon));
    
    
    yRescaled = normalize(xlat(:), 1, 'range', [1, size(xlon,2)]);
    yRescaled = reshape(yRescaled, size(xlon));
    
    d = cat(3, (yRescaled-py), (xRescaled-px));
     
    im = imwarp(var,-d);
    
    alphaData = ones(size(im))*0.8;
    alphaData(im == 0) = 0;
    alphaData(isnan(im)) = 0;
    
    % close all
    imsc=imagesc(ax,...
        [min(xlon(:)), max(xlon(:))], ...
        [min(xlat(:)), max(xlat(:))], ...
        im', ...
        "AlphaData", alphaData' ...
    );
    set(ax, 'YDir','normal') 
end

