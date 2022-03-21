function [ C_refined, matches_refined ] = refine( Src, Tar, C, varargin )
%REFINE Refines the Fmap
    p = inputParser;
    addParameter(p,'type','ICP',@ischar);
    default_params = struct;
    default_params.numIter = 5;
    default_params.stepSize = 1;
    default_params.startSize = 10;
    addParameter(p,'parameters',default_params,@isstruct);
    parse(p,varargin{:});
    inputs = p.Results;
    
    type = inputs.type;
    parameters = inputs.parameters;
    
    switch type
        case 'ICP'
            assert(isfield(parameters, 'numIter'),...
                'parameters have to contain field numIter when using ICP refinement.');
            
            if isfield(Src, 'BASIS_ICP') && isfield(Tar, 'BASIS_ICP')
                BASIS_Src = Src.BASIS_ICP;
                BASIS_Tar = Tar.BASIS_ICP;
                % Converting C to new basis
                C = BASIS_Tar.basis_inverse*Tar.BASIS.basis*C*Src.BASIS.basis_inverse*BASIS_Src.basis;
            else
                BASIS_Src = Src.BASIS;
                BASIS_Tar = Tar.BASIS;
            end
            
            [C_refined, matches_refined] = icp_refine(BASIS_Src, BASIS_Tar, C, parameters.numIter);  
            % Converting C back to the original basis
            C_refined = Tar.BASIS.basis_inverse*BASIS_Tar.basis*C_refined*BASIS_Src.basis_inverse*Src.BASIS.basis;
        case 'ZoomOut'
            assert(isfield(parameters, 'stepSize'),...
                'parameters have to contain field stepSize when using ZoomOut refinement.');
            assert(isfield(parameters, 'startSize'),...
                'parameters have to contain field startSize when using ZoomOut refinement.');
            
            matches = knnsearch((C*Src.BASIS.basis')',Tar.BASIS.basis);
            
            if isfield(Src, 'BASIS_ZoomOut') && isfield(Tar, 'BASIS_ZoomOut')
                BASIS_Src = Src.BASIS_ZoomOut;
                BASIS_Tar = Tar.BASIS_ZoomOut;
            else
                BASIS_Src = Src.BASIS;
                BASIS_Tar = Tar.BASIS;
            end
            
            matches_refined = zoomOut_refine(BASIS_Src, BASIS_Tar, ...
                matches, ...
                parameters.startSize,parameters.stepSize);
            
            C_refined = BASIS_Tar.basis\BASIS_Src.basis(matches_refined,:);
    end
end

