function [ NN ] = knn( tar, src )
    %KNN Wrapper around annquery to reproduce knnsearch from the Statistics
    %and Machine Learning Toolbox of Matlab.
    
    NN = double(transpose(annquery(transpose(tar),transpose(src),1)));
end

