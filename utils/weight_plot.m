function [ ] = weight_plot( A, X, colormap, edge_width )

    default('colormap', winter);
    default('edge_width', 2); 
    %Do log weights to get a better spread for the color map.
    A_log = log(A) + .001;
    A_log(A_log == -inf) = 0;

    wgPlot(A_log,X,'edgeColorMap',colormap, 'edgewidth', edge_width);

end

