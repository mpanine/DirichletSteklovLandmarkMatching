function [fullp2pTarSrc, fullp2pSrcTar, FmapSrcTar_combined, FmapTarSrc_combined] = findFullP2P_Iterative(Src_eF, Src_evals, SrcMass, Tar_eF, Tar_evals, TarMass, FmapSrcTar, FmapTarSrc)
%Converts the functional maps obtained at the boundary into point to point
%maps in the bulk. Does so iteratively, like a ZoomOut with no increase in
%basis size.

%% Internal settings

toplot = true; % Set to true for figures of the Fmaps.
numit = 10;

%%


num_bounds = size(Src_eF,3);
num_eigs = size(Src_eF,2);


%% Pre-allocation
Src_evecs_combined = zeros(size(Src_eF,1), size(Src_eF,2)*size(Src_eF,3));
Tar_evecs_combined = zeros(size(Tar_eF,1), size(Tar_eF,2)*size(Tar_eF,3));

Src_evals_combined = zeros(num_bounds*num_eigs, num_bounds*num_eigs);
Tar_evals_combined = zeros(num_bounds*num_eigs, num_bounds*num_eigs);

FmapSrcTar_combined = zeros(num_bounds*num_eigs, num_bounds*num_eigs);
FmapTarSrc_combined = zeros(num_bounds*num_eigs, num_bounds*num_eigs);
% TEST THIS IF ABOVE IS SLOW
% Src_evals_combined = spalloc(num_bounds*num_eigs, num_bounds*num_eigs, num_bounds*num_eigs);
% Tar_evals_combined = spalloc(num_bounds*num_eigs, num_bounds*num_eigs, num_bounds*num_eigs);
% 
% FmapSrcTar_combined = spalloc(num_bounds*num_eigs, num_bounds*num_eigs, num_bounds*num_eigs*num_bounds*num_eigs);
% FmapTarSrc_combined = spalloc(num_bounds*num_eigs, num_bounds*num_eigs, num_bounds*num_eigs*num_bounds*num_eigs);


for i = 1:num_bounds % Combining the arrays
        
       start_ind = (i-1)*num_eigs + 1;
       end_ind = i*num_eigs;
       

       
       Src_evecs_combined(:, start_ind:end_ind ) = Src_eF(:,:,i);
       Tar_evecs_combined(:, start_ind:end_ind ) = Tar_eF(:,:,i);
       

       
       Src_evals_combined(start_ind:end_ind, start_ind:end_ind) = Src_evals(:,:,i);
       Tar_evals_combined(start_ind:end_ind, start_ind:end_ind) = Tar_evals(:,:,i);
    
       FmapSrcTar_combined(start_ind:end_ind, start_ind:end_ind) = FmapSrcTar{i};
       FmapTarSrc_combined(start_ind:end_ind, start_ind:end_ind) = FmapTarSrc{i};
end


%% Pseudo-inverting the comobined Dirichlet-Steklov basis (DOES THIS MAKE SENSE???)

TarProjOnDS = pinv(Tar_evecs_combined);
SrcProjOnDS = pinv(Src_evecs_combined);

% figure
% imagesc(abs(TarProjOnDS*Tar_evecs_combined))
% figure
% imagesc(abs(SrcProjOnDS*Src_evecs_combined))


%% Plot the original Fmaps

if toplot
    
    figure('Units','normalized','Position',[0.1 0.1 0.8 0.8])
    subplot(1,3,1)
    imagesc(FmapSrcTar_combined); title('Fmap Src->Tar; It: 0')
%     axis equal
    
    subplot(1,3,2)
    imagesc(FmapTarSrc_combined); title('Fmap Tar->Src; It: 0')
%     axis equal
    
    subplot(1,3,3)
    imagesc(FmapSrcTar_combined*FmapTarSrc_combined); title('Bijectivity check; It: 0')
%     axis equal
    
    pause(0.1)
    
end


for i = 1:numit



    %% Computing the fullp2pTarSrc map.

        SrcFTarSrc = Src_evecs_combined*Src_evals_combined +... %Steklov term
               Src_evecs_combined+... %"Proper" functional map term (pullback of p2p map)
               Src_evecs_combined*FmapSrcTar_combined'+... %Orthonormality term.
               Src_evecs_combined*FmapTarSrc_combined; %Bijectivity term

        TarFSrcTar = Tar_evecs_combined*Tar_evals_combined*FmapSrcTar_combined +... %Steklov term
               Tar_evecs_combined*FmapSrcTar_combined+... %"Proper" functional map term (pullback of p2p map)
               Tar_evecs_combined+...%Orthonormality term.
               Tar_evecs_combined; %Bijectivity term


        fullp2pTarSrc = annquery(SrcFTarSrc', TarFSrcTar', 1);



        %% Computing the fullp2pSrcTar map.

        SrcFSrcTar = Src_evecs_combined*Src_evals_combined*FmapTarSrc_combined +... %Steklov term
               Src_evecs_combined*FmapTarSrc_combined+... %"Proper" functional map term (pullback of p2p map)
               Src_evecs_combined+... %Orthonormality term.
               Src_evecs_combined; %Bijectivity term

        TarFTarSrc = Tar_evecs_combined*Tar_evals_combined +... %Steklov term
               Tar_evecs_combined+... %"Proper" functional map term (pullback of p2p map)
               Tar_evecs_combined*FmapTarSrc_combined'+...%Orthonormality term.
               Tar_evecs_combined*FmapSrcTar_combined; %Bijectivity term


        fullp2pSrcTar = annquery(TarFTarSrc', SrcFSrcTar', 1);


        %% Compute the new Fmaps
        
%         FmapSrcTar_combined = Tar_evecs_combined' * TarMass * Src_evecs_combined(fullp2pTarSrc,:); %WRONG INNER PRODUCT
%         FmapTarSrc_combined = Src_evecs_combined' * SrcMass * Tar_evecs_combined(fullp2pSrcTar,:);
        
        FmapSrcTar_combined = TarProjOnDS * Src_evecs_combined(fullp2pTarSrc,:);        
        FmapTarSrc_combined = SrcProjOnDS * Tar_evecs_combined(fullp2pSrcTar,:);
        
        
        %% Plot the new Fmaps
        
        if toplot
            
            clf
            
            subplot(1,3,1)
            imagesc(abs(FmapSrcTar_combined)); title(sprintf('Fmap Src->Tar; It: %d',i))
%             axis equal

            subplot(1,3,2)
            imagesc(abs(FmapTarSrc_combined)); title(sprintf('Fmap Tar->Src; It: %d',i))
%             axis equal

            subplot(1,3,3)
            imagesc(abs(FmapSrcTar_combined*FmapTarSrc_combined)); title(sprintf('Bijectivity check; It: %d',i))
%             axis equal

            pause(1e-9)

        end
        

end
















end