function writeMesh_plyModify(fileName, vertices, faces, color, rgb)
% Optimized WRITEMESH_PLY to write a mesh into a text file in PLY format.

%% Check inputs
if ~ischar(fileName)
    error('First argument must contain the name of the file');
end

% Optionally parse struct data
if isstruct(vertices)
    faces = vertices.faces;
    vertices = vertices.vertices;
end

%% Initializations
cmin = min(color(color~=0));
cmax = max(color);

nVertices = size(vertices, 1);
nFaces = size(faces, 1);
if iscell(faces)
    nFaces = length(faces);
end

% Open file for writing text
f = fopen(fileName, 'wt');
if (f == -1)
    error('Couldn''t open the file %s', fileName);
end

% Prepare colormap and color scaling if `rgb` is not provided
if nargin < 5
    colorMap = 'jet';
    numBins = 512;
    
    % Use built-in colormap generation
    cMap = colormap(strcat(colorMap, '(', num2str(numBins), ')'));
    
    % Discretize the color data and map to colormap
    binEdges = linspace(cmin, cmax, numBins + 1);
    binEdges(1) = min(min(color), binEdges(1));
    binEdges(end) = max(max(color), binEdges(end));
    [cVectorDiscreet, ~] = discretize(color, binEdges);
    
    % Map discrete color values to RGB format (0-255 range)
    c = cMap(cVectorDiscreet, :); 
    colors = round(rescale(c, 0, 255));
else
    colors = round(rescale(rgb, 0, 255));
end

%% Write Header
fprintf(f, 'ply\n');
fprintf(f, 'format ascii 1.0\n');
fprintf(f, 'comment created by MatGeom for Matlab\n');
fprintf(f, 'element vertex %d\n', nVertices);
fprintf(f, 'property float x\nproperty float y\nproperty float z\n');
fprintf(f, 'element face %d\n', nFaces);
fprintf(f, 'property list uchar int vertex_index\n');
fprintf(f, 'property uchar red\nproperty uchar green\nproperty uchar blue\n');
fprintf(f, 'end_header\n');

%% Write Vertex Data (Vectorized)
% Write vertices in one go
fprintf(f, '%g %g %g\n', vertices');

%% Write Face Data (Vectorized)
% Preconstruct the face and color data for efficient writing
faceData = [repmat(3, nFaces, 1), faces - 1, colors];

% Write faces and their colors in one go
fprintf(f, '%d %d %d %d %d %d %d\n', faceData');

% Close file
fclose(f);
end

