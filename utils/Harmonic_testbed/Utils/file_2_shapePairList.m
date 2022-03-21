function [ shapePairList ] = file_2_shapePairList( filename )
    %FILE_2_SHAPEPAIRLIST 
    fid = fopen(filename,'rt');
    src = textscan(fid,'%s','delimiter','\n');
    fclose(fid);
    shapePairList = src{1}';
end

