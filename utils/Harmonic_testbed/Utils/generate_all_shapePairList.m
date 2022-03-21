function [ shapePairList ] = generate_all_shapePairList( shapeList )
    %GENERATE_ALL_SHAPEPAIRLIST Generates all shape pairs given a list of
    %shapes.
    
    shapePairList = repmat({''},1,size(shapeList, 2)*size(shapeList, 2));
    for i=1:size(shapeList, 2)
        for j=1:size(shapeList, 2)
            shapePairList{i+(j-1)*size(shapeList, 2)} = [shapeList{j},' ',shapeList{i}];
        end
    end
    
end

