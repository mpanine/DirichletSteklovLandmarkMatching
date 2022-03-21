function [ outputs ] = compute_pMap_evaluations( Src, Tar, ...
                                            Src_eval_landmarks, Tar_eval_landmarks, ...
                                            pF_SrcTar, varargin )
    %COMPUTE_PMAP_EVALUATIONS Computes the evaluation of pF_SrcTar based on
    %the evaluation methods specified in varargin:
    % methods:
    % - geodesic error: 'geodesic'
    % - euclidean error: 'euclidean'
    % - distortion: 'distortion'
    % - dirichlet: 'dirichlet'
    % geodesic error by default:
    if nargin<6
        varargin = {'geodesic'};
    end;
    
    outputs = {};
    for i=1:length(varargin)
        cur = char(varargin{i});
        fprintf(sprintf('%s...',cur));
        
        switch cur
            case 'geodesic'
                outputs{end+1} = evaluation.geodesic_metrics(Tar, pF_SrcTar, Src_eval_landmarks, Tar_eval_landmarks);
                fprintf(sprintf('Mean error = %d',outputs{end}.mean));
            case 'euclidean'
                outputs{end+1} = evaluation.euclidean_metrics(Tar, pF_SrcTar, Src_eval_landmarks, Tar_eval_landmarks);
            case 'distortion'
                outputs{end+1} = evaluation.distortion_metrics(Src, Tar, pF_SrcTar); 
            case 'dirichlet'
                outputs{end+1} = evaluation.map_dirichlet_energy(Src, Tar, pF_SrcTar); 
        end;
    end;
    
end

