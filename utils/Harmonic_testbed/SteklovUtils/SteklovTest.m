% Temporary test for the Steklov code.

clear all; close all;

shape_settings = shape.Settings('computeMHB', true,...
                        'computeGeoDist', false,...
                        'MHB_settings', basis.MHB_Settings(20),...
                        'computeNormals', false);
cache_settings = cache.Settings('a','b', false);

Src = struct;
Src.SHAPE = shape.compute('./Data/Faust_remeshed/vtx_5k/tr_reg_000.off', ...
                        shape_settings, cache_settings);
                    
% Src.SHAPE = shape.compute('./Data/bunny_conformal/flow00000000.off', ...
%                         shape_settings, cache_settings);


centerpt = 4000;% Positions the segment %Faust 4999 + thresh = 0.25 : lots of boundary components
thresh = 0.25;%0.25;%0.3;% Size of the segment (eucliean dist from centerpt) 0.49 on Faust: touching boundary components

ifDirichletAtCenter = true; %Sets a homogeneous Dirichlet BC at the center if true.



dists = sqrt( sum( (Src.SHAPE.surface.VERT - Src.SHAPE.surface.VERT(centerpt,:)).^2 , 2));

segment = find(dists <= thresh);

% [boundary_old, boundary_new, boundary_edges_old] = findSegmentBoundary(Src.SHAPE, segment);
[W, S, segment_TRIV, boundary_old, boundary_new, boundary_edges_old] = splitSegmentSteklov(Src.SHAPE, segment);

new_center = find(segment == centerpt);

WW = W; %This is not especially elegant;
SS = S;

if ifDirichletAtCenter
    
    
    
    WW(new_center,:) = [];
    WW(:,new_center) = [];
    
    SS(new_center,:) = [];
    SS(:,new_center) = [];
    
end

                  
num_eigs = min(25,length(boundary_new));

while num_eigs > 0 %This loop reduces the number of requested eigenfunctions until we get real eigenvectors.

    try
        [evecs, evals] = eigs(WW, SS, num_eigs, 1e-5);
    catch
        % In case of trouble make the laplacian definite
        [evecs, evals] = eigs(WW - 1e-8*speye(size(WW)), SS, num_eigs, 'sm');
        warning('The Steklov boundary is likely bad.')
    end
    
	prods = evecs'*SS*evecs;
    inv_pr = diag( diag(prods).^-0.5 ); %Normalization, just in case (sometimes necessary) (WHY?)
    evecs = evecs * inv_pr;
    
    if isreal(evecs)
        break;
    else
        num_eigs = num_eigs - 5;
    end
    
end

if ifDirichletAtCenter
    
    evecs = reinsertRemovedPoint(evecs, new_center);
    
end



eF = extendSegmentFunctions(evecs, segment, Src.SHAPE.nv);

%% Prepare boundary for plotting
boundary_reorder = boundaryReorder(boundary_edges_old);
reorder_lengths = computeLengthAlongBoundary(boundary_reorder, Src.SHAPE);

%% Plot the overall shape

plotShapeAndLandmarks(Src.SHAPE, segment, boundary_old)
title('Shape, segment and boundary')
hold on

plot3(Src.SHAPE.surface.X(centerpt),Src.SHAPE.surface.Y(centerpt),Src.SHAPE.surface.Z(centerpt),'.','Color','k','MarkerSize',20)

for i = 1:size(boundary_edges_old,1) % Plot all edges separately: no need to worry about connectivity
    
    X = Src.SHAPE.surface.VERT(boundary_edges_old(i,:),1);
    Y = Src.SHAPE.surface.VERT(boundary_edges_old(i,:),2);
    Z = Src.SHAPE.surface.VERT(boundary_edges_old(i,:),3);
    
    plot3(X,Y,Z,'k')
    
end



% for i = 1:length(boundary_reorder)
%     
%     % Must close loops if appropriate
%     %TODO: make sure this covers all possibili
%     edge = [boundary_reorder{i}(1) boundary_reorder{i}(end)];    
%     edge_ind = find(ismember(boundary_edges_old, edge, 'rows'));
%     if isempty(edge_ind)
%         edge_ind = find(ismember(boundary_edges_old, circshift(edge,1),'rows'));
%     end
%     
%     if ~isempty(edge_ind)
%         bb = [boundary_reorder{i} boundary_reorder{i}(1)]; %Closed loop
%     else
%         bb = boundary_reorder{i}; %Open ends
%     end
%     
%     
%     
%     
%     X = Src.SHAPE.surface.VERT(bb,1);
%     Y = Src.SHAPE.surface.VERT(bb,2);
%     Z = Src.SHAPE.surface.VERT(bb,3);
%     
%     plot3(X,Y,Z,'k')
%     
% end


%% Plot the segment by itself

% figure
% trimesh(segment_TRIV, Src.SHAPE.surface.X, Src.SHAPE.surface.Y, Src.SHAPE.surface.Z, ...
% 'FaceColor','interp', ...
% 'FaceAlpha', 1, 'EdgeColor', [0.9 0.9 0.9]); axis equal;
% title('The segment by itself')


%% Spectrum and Normalization check

figure
subplot(1,2,1)
plot(diag(evals),'.')
title('Steklov Eigenvalues')
subplot(1,2,2)
imagesc(evecs'*S*evecs); colorbar;
title('Orthogonalization check')

%% Plot boundary eigenfunctions

funs_to_plot = 1:2;



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

%%
nef = 20;
figure
subplot(1,2,1)
trimesh(segment_TRIV, Src.SHAPE.surface.X, Src.SHAPE.surface.Y, Src.SHAPE.surface.Z, eF(:,nef), ...
'FaceColor','interp', ...
'FaceAlpha', 1, 'EdgeColor', 'none'); axis equal; axis off; colorbar;
title(sprintf('Steklov Eigenfunction %d',nef))

subplot(1,2,2)
trimesh(segment_TRIV, Src.SHAPE.surface.X, Src.SHAPE.surface.Y, Src.SHAPE.surface.Z, Src.SHAPE.W*eF(:,nef), ...
'FaceColor','interp', ...
'FaceAlpha', 1, 'EdgeColor', 'none'); axis equal; axis off; colorbar;
title(sprintf('Laplacian of Same',nef))



