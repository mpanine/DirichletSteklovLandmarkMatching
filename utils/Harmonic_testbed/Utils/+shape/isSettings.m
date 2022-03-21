function [ is_class ] = isSettings( input_obj )
%ISSETTINGS Returns whether input is an shape.Settings object
    is_class = isa(input_obj, 'shape.Settings');
end


