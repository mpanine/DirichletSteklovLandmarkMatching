function [  ] = save_pMap_evaluations( evaluation_outputs, method, i, Src_name, Tar_name, saveDir, varargin )
    %SAVE_PMAP_EVALUATIONS 
    % /!\ varargin order has to be consistent with varargin of
    % compute_pMap_evaluations /!\
    if nargin<6
        varargin = {'geodesic'};
    end;
    
    assert(length(varargin)==length(evaluation_outputs));
    
    for j=1:length(varargin)
        cur = char(varargin{j});
        saveName = sprintf('%s/%d_%s-%s_%s.mat',char(saveDir), i, Src_name, Tar_name,cur);
        
        computed_evals = struct;
        if exist(saveName,'file')
            computed_evals = load(saveName);
        end;
        
        if isfield(computed_evals,[method,'_',cur])
            warning('Shape pair %s-%s already evaluated for method %s',shape_name_Src,shape_name_Tar,cur_method);
        else
            result.([method,'_',cur]) = evaluation_outputs{j};

            if exist(saveName, 'file')
                save(saveName,'-struct','result','-append');
            else
                save(saveName,'-struct','result');
            end;
        end;
    end;
    
end

