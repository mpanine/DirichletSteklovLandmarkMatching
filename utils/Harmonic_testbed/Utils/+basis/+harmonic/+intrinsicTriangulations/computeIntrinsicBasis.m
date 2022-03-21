function Basis = computeIntrinsicBasis( shape_number, numvertices, landmarks, landmarks_radii )
%computeIntrinsicBasis computes the intrinsic basis restricted to the input vertices
%shape_number: number of the shape. This refers to the global variable storing the path to the .obj version of the shape.
%numvertices: number of vertices in the initial shape
%landmarks: indices of the landmarks (NOT INCLUDING CONTROL LANDMARKS)
%landmarks_radii: list of the radii of the circumcircles of the triangles near adjacent to the landmark


minrad = 0.000001; % Smallest allowed landmark radius. Useful in an optimization code which could output nonpozitive radii.
if min(landmarks_radii) < minrad
    
    warning('The at least one landmark radius is below minrad. It has been set to minrad.')
    landmarks_radii( landmarks_radii < minrad) = minrad;   
    
end


global objpath1;
global objpath2;

if shape_number == 1
    path = objpath1;
elseif shape_number == 2
    path = objpath2;
else
    error('The shape number is not either 1 or 2.')
end

landmarks = landmarks(:)';
landmarks_radii = landmarks_radii(:)';


prefile = [landmarks - 1; landmarks_radii]; % -1 required for 0-indexing

fileID = fopen('refineAroundVerts.txt','w');
fprintf(fileID, '%d %f\n', prefile);
fclose(fileID);

%! ./IntrinsicTriangulations/navigating-intrinsic-triangulations-demo/build/bin/int_tri  path
%--noGUI --refineAroundVerts=[refineAroundVerts.txt] --laplaceMat
system(['./IntrinsicTriangulations/navigating-intrinsic-triangulations-demo/build/bin/int_tri '  path ' --noGUI --refineAroundVerts ./refineAroundVerts.txt --laplaceMat --interpolateMat']);


fid = fopen('./laplace.spmat') ;
preW = textscan(fid,'%f %f %f','HeaderLines',1) ;
preW = cell2mat(preW);
fclose(fid);

W = spconvert(preW);

% imagesc(W)
% size(W)

preBasis = ComputeHarmonicBasisIntrinsic(W, landmarks);

Basis = preBasis(1:numvertices,:);

% fid = fopen('./interpolate.spmat') ;
% preInt = textscan(fid,'%f %f %f','HeaderLines',1) ;
% preInt = cell2mat(preInt);
% fclose(fid);
% 
% Int = spconvert(preInt)











end