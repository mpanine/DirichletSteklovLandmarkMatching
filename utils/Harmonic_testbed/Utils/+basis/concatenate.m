function [ BASIS_CONCAT ] = concatenate( BASIS_ls )
    %CONCATENATE Concatenate the bases contained in BASIS_ls in the order
    %of BASIS_ls.
    
    BASIS_CONCAT = struct;
    BASIS_CONCAT.basis = BASIS_ls{1}.basis;
    BASIS_CONCAT.basis_inverse = BASIS_ls{1}.basis_inverse;
    BASIS_CONCAT.settings = {BASIS_ls{1}.settings};
    for i=2:length(BASIS_ls)
        BASIS_CONCAT.basis = [BASIS_CONCAT.basis,BASIS_ls{i}.basis];
        BASIS_CONCAT.basis_inverse = [BASIS_CONCAT.basis_inverse;BASIS_ls{i}.basis_inverse];
        BASIS_CONCAT.settings{end+1} = BASIS_ls{i}.settings;
    end
end

