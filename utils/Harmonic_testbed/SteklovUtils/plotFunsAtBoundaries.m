function plotFunsAtBoundaries(eF, boundary_reorder, reorder_lengths, funs_to_plot)
%Plots eF at the reordered Steklov boundary
% the functions to plot are indexed in funs_to_plot


num_boundaries = length(boundary_reorder);

ymax = -Inf; ymin = Inf;

for i = 1:num_boundaries % Find the maxima and minima of the function to plot.
    
    yma = max(max( eF(boundary_reorder{i},funs_to_plot) ));
    ymi = min(min( eF(boundary_reorder{i},funs_to_plot) ));
    
    if yma > ymax
        ymax = yma;
    end
    
    if ymi < ymin
       ymin = ymi; 
    end

end

yy = max(abs(ymax),abs(ymin));

ymax = ymax + 0.05*abs(yy);
ymin = ymin - 0.05*abs(yy);

figure
for i = 1:num_boundaries
    
    subplot(1, num_boundaries, i)
    plot(reorder_lengths{i}, eF(boundary_reorder{i},funs_to_plot) )
    if num_boundaries > 1
        title(sprintf('Component %d',i))
    end
    xlim([min(reorder_lengths{i}) max(reorder_lengths{i})])
    ylim([ymin ymax])
    xlabel('Arclength along boundary')

end

for i = 1:length(funs_to_plot)
   
    ll{i} = sprintf('Eigenfunction %d',funs_to_plot(i));
    
end
legend(ll)
% sgtitle('Steklov Eigenfunctions on the Steklov boundary')





end