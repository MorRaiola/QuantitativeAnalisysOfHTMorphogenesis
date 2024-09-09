function writeMesh_plyModify(fileName, vertices, faces,color,rgb)
%WRITEMESH_PLY Writes a mesh into a text file in PLY format.
%
%   writeMesh_ply(FNAME, VERTS, FACES)
%
%   Example
%   writeMesh_ply
%
%   See also
%   meshes3d, readMesh_ply, writeMesh_off
 
% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2018-04-26,    using Matlab 9.4.0.813654 (R2018a)
% Copyright 2018 INRA - Cepia Software Platform.


%% Check inputs

if ~ischar(fileName)
    error('First argument must contain the name of the file');
end

% optionnaly parses data
if isstruct(vertices)
    faces = vertices.faces;
    vertices = vertices.vertices;
end


cmin = min(color(color~=0));
cmax = max (color);

%% Initializations

% number of vertices and faces
nVertices = size(vertices, 1);
nFaces = size(faces, 1);
if iscell(faces)
    nFaces = length(faces);
end

% open file for writing text
f = fopen(fileName, 'wt');
if (f == -1)
	error('Couldn''t open the file %s', fileName);
end


%% Write Header 

% write the header line
fprintf(f, 'ply\n');

% write format (only ASCII supported)
fprintf(f, 'format ascii 1.0\n');

% some comments
fprintf(f, 'comment created by MatGeom for Matlab\n');


% valori = unique(color);
% 
% RGB = jet(size(valori,1));
% colors= zeros(size(color,1),3);
% for h=2:size(valori,1)
%     colors(color ==valori(h),:)= repmat(RGB(h,:),size(colors(color ==valori(h),1)));
% end
% 
% colors = round(rescale(colors,0,255));
% % colors = 255*ones(nFaces,3);

if nargin<5
colorMap = ('jet'); RGB_bins = 512;
cMap=colormap(strcat(colorMap, '(', num2str(RGB_bins), ')'));
cbin_edges=linspace(cmin, cmax, RGB_bins+1);
    if min(color)<cbin_edges(1)
        cbin_edges(1)=min(color);  
    end
    if max(color)>cbin_edges(end)
        cbin_edges(end)=max(color);
    end
[cVector_discreet, ~]=discretize(color, cbin_edges);
c = cMap(cVector_discreet, :); 
colors = round(rescale(c,0,255));
else
colors = round(rescale(rgb,0,255));
end

% write declaration for vertices
fprintf(f, 'element vertex %d\n', nVertices);
fprintf(f, 'property float x\nproperty float y\nproperty float z\n');


% write declaration for faces
fprintf(f, 'element face %d\n', nFaces);
fprintf(f, 'property list uchar int vertex_index\n');
fprintf(f, 'property uchar red\nproperty uchar green\nproperty uchar blue\n');

% end of header
fprintf(f, 'end_header\n');

%% Write vertex info

format = '%g %g %g\n';
for iv = 1:nVertices
    fprintf(f, format, vertices(iv, :));
end


%% Write face info
for iFace = 1:nFaces
    fprintf(f, '3 %d %d %d ', faces(iFace, :) - 1);
    fprintf(f, '%g %g %g\n', colors(iFace,:));
end
% close the file
fclose(f);
close
end
