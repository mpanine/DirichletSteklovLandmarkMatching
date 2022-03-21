function [ grad_colors ] = color_gradient(dark, length)
    white = [1, 1, 1];
    grad_colors = [linspace(dark(1),white(1),length)',...
              linspace(dark(2),white(2),length)',...
              linspace(dark(3),white(3),length)'];
end

