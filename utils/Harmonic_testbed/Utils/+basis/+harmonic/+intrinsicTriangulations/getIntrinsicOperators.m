function [W, A] = getIntrinsicOperators( shape_path, landmarks, landmarks_radii )
%getIntrinsicOperators outputs the stiffness and mass matrices after an intrinsic refinement
%shape_path: path of the .obj version of the shape.
%landmarks: indices of the landmarks (NOT INCLUDING CONTROL LANDMARKS, if any)
%landmarks_radii: list of the radii of the circumcircles of the triangles near adjacent to the landmark
import basis.harmonic.intrinsicTriangulations.*

warning('error', 'MATLAB:DELETE:FileNotFound');
warning('error', 'MATLAB:DELETE:Permission');
fprintf('\n');
try
    delete 'refineAroundVerts.txt' 
catch
end
try
    delete 'faceInds.dmat' 
catch
end
try 
    delete 'vertexPositions.dmat' 
catch
end
try 
    delete 'laplace.spmat' 
catch
end
try 
    delete 'faceLengths.dmat'
catch
end

    

% minrad = 0.000001; % Smallest allowed landmark radius. Useful in an optimization code which could output nonpositive radii.
% if min(landmarks_radii) < minrad
%     
%     warning('At least one landmark radius is below minrad. It has been set to minrad.')
%     landmarks_radii( landmarks_radii < minrad) = minrad;   
%     
% end


landmarks = landmarks(:)';
landmarks_radii = landmarks_radii(:)';


prefile = [landmarks - 1; landmarks_radii]; % -1 required for 0-indexing

fileID = fopen([pwd '\refineAroundVerts.txt' ],'w');
fprintf(fileID, '%d %f\n', prefile);
fclose(fileID);

%! ./IntrinsicTriangulations/navigating-intrinsic-triangulations-demo/build/bin/int_tri  path
%--noGUI --refineAroundVerts=[refineAroundVerts.txt] --laplaceMat

%Linux Path
% system(['./harmonic/IntrinsicTriangulations/navigating-intrinsic-triangulations-demo/build/bin/int_tri '  shape_path ' --noGUI --refineDelaunayAroundVerts ./refineAroundVerts.txt --intrinsicFaces --vertexPositions --laplaceMat']);

%Windows Path

system([pwd '\utils\+basis\+harmonic\+intrinsicTriangulations\navigating-intrinsic-triangulations-demo\build\bin\int_tri.exe '  pwd '\' shape_path ' --noGUI --refineDelaunayAroundVerts \refineAroundVerts.txt --intrinsicFaces --vertexPositions --laplaceMat']);


fid = fopen([pwd '\faceInds.dmat']) ;
preTRIV = textscan(fid,'%d %d %d','HeaderLines',1) ;
TRIV = cell2mat(preTRIV);
TRIV = TRIV + 1;% 1-index instead of 0-index.
fclose(fid);

fid = fopen([pwd '\vertexPositions.dmat']) ;
preX = textscan(fid,'%f %f %f','HeaderLines',1) ;
X = cell2mat(preX);
fclose(fid);


fid = fopen([pwd '\faceLengths.dmat']) ;
preEdgeLengths = textscan(fid,'%f %f %f','HeaderLines',1) ;
EdgeLengths = cell2mat(preEdgeLengths);
fclose(fid);


% [W, A] = cotLaplacian_for_Intrinsic(X, TRIV);
[~, A] = cotLaplacian_for_Intrinsic(EdgeLengths, TRIV);

% A = sparse(1:length(A), 1:length(A), A);
area = sum(A);
A = sparse(diag(A)/area);

fid = fopen([pwd '\laplace.spmat']) ;
preW = textscan(fid,'%f %f %f','HeaderLines',1) ;
preW = cell2mat(preW);
fclose(fid);

W = spconvert(preW);



end