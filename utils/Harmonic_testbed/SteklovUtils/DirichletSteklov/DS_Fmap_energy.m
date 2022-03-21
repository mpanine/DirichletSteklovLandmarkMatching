function E = DS_Fmap_energy(Src_evals, Tar_evals, Fmap_SrcTar, Fmap_TarSrc)
%DS_Fmap_energy: computes the energy of a Dirichlet-Steklov Fmap: Src -> Tar
%Used in the optimization of the initial correspondence.
%Works on one boundary component.

E = norm(  Fmap_SrcTar  * Src_evals    -   Tar_evals * Fmap_SrcTar   )^2 + ... %Commutativity with Steklov
    norm(  Fmap_SrcTar' * Fmap_SrcTar  -   eye(size(Fmap_SrcTar))    )^2 + ... %Othonormality part
    norm(  Fmap_TarSrc  * Fmap_SrcTar  -   eye(size(Fmap_SrcTar))    )^2 ;     %Bijectivity part


end