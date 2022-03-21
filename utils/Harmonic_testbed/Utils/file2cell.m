function [ results ] = file2cell( directory, name, entry_type )
    %FILE2CELL Converts folder containing stored results to a cell array.
    
    saveDir = [directory,name];
    
    set = load([saveDir,'/settings.mat']);
    
    if ~exist(saveDir,'dir')
        throw(MException('FILE2CELLERROR:directoryUnknown',...
                         'The directory directory/name/ does not exist.'));
    end
    
    if exist('set.singleBasisEvaluation_settings','var')
        s = set.singleBasisEvaluation_settings;
    elseif isfield(set,'shapePairList')
        s = set;
    else
        throw(MException('FILE2CELLERROR:noShapePairList',...
                         'The file directory/name/settings.mat does not contain a shapePairList variable.'));
    end
    
    entries = {};
    
    for i=1:length(s.shapePairList)
        fprintf(sprintf('Iteration %d/%d\n',i,length(s.shapePairList)));
        pair_str = s.shapePairList(i);
        [names, ~] = split(pair_str);
        shape_name_Src = char(names(1));
        shape_name_Tar = char(names(2));
        
        try
            entries{end+1} = load(sprintf('%s/%d_%s-%s_%s.mat',char(saveDir), i, shape_name_Src, shape_name_Tar, entry_type));
        catch
            warning(sprintf('Problem loading the data at iteration %d (shape pair %s-%s)',i,shape_name_Src,shape_name_Tar));
            break
        end
        
        clear result
    end;
    
    results.entries = entries;
    
end

