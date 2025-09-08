%% Validating Spatial Correspondences between ATLAS and Live-Shape

% Define paths and filenames
liveShapeFile = '\Code_data\2. IntegratingMultipleLiveImagesIntoAConsensusTemporalReference\Cluster5\3D\Embryo\27.ply';
atlasFile = '\Code_data\2. IntegratingMultipleLiveImagesIntoAConsensusTemporalReference\Cluster5\3D\Atlas\27.ply';
correspondenceFile = '\Code_data\2. IntegratingMultipleLiveImagesIntoAConsensusTemporalReference\Cluster5\VertexCorrespondence\Embryo27\IdxCUT.mat';
matchingFile = '\Code_data\2. IntegratingMultipleLiveImagesIntoAConsensusTemporalReference\Cluster5\VertexCorrespondence\Embryo27\IdxMatch.mat';
savePath = '\Code_data\2. IntegratingMultipleLiveImagesIntoAConsensusTemporalReference\Cluster5\Deformation\';
outputFile = 'indexSmooth27.mat';

% Load live shape data
[node, face] = read_ply(liveShapeFile);
centroid = meshcentroid(node, face);

% Initialize geodesic library and settings
global geodesic_library;                
geodesic_library = 'geodesic_release'; 


% Initialize new mesh and geodesic algorithm
mesh = geodesic_new_mesh(node, face); 
algorithm = geodesic_new_algorithm(mesh, 'exact'); 

% Preallocate index cell array to store valid face indices
index = cell(size(face, 1), 1);

% Loop through each face of the mesh
for i = 1:size(face, 1)
    % Find the nearest vertex to the current face centroid
    vertex_id = knnsearch(node, centroid(i, :)); 

    % Create a single source at the identified vertex
    source_points = {geodesic_create_surface_point('vertex', vertex_id, node(vertex_id, :))};

    % Propagate the geodesic algorithm from the source point
    geodesic_propagate(algorithm, source_points, [], 20); 

    % Retrieve distances to all vertices
    [~, distances] = geodesic_distance_and_source(algorithm); 

    % Identify vertices within a valid distance range
    valid_vertices = find(distances < 1e10);
    
    % Initialize a logical array to mark valid faces
    valid_face = false(size(face, 1), 1);
    
    % Check if all vertices of each face are within the valid distance
    for d = 1:size(face, 1)
        if all(ismember(face(d, :), valid_vertices))
            valid_face(d) = true;
        end
    end
    
    % Store the indices of valid faces
    index{i} = valid_face;
end


% Save the index data for smooth regions
save(fullfile(savePath, outputFile), 'index');

%% Smooth the Area and Face-to-Face Matching

% Load ATLAS shape data
[nodeA, faceA] = read_ply(atlasFile);

% Load correspondence and matching data
load(correspondenceFile, 'IdxCUT');
load(matchingFile, 'IdxMatch');

% Load the saved index file
load(fullfile(savePath, outputFile), 'index');

% Calculate face areas for the live-shape
area = meshFaceAreas(node, face);

% Calculate the smoothed area values based on valid faces
newarea = zeros(size(face, 1), 1);
for z = 1:size(face, 1)
    pos = find(index{z});
    if ~isempty(pos)
        newarea(z) = mean(area(pos));
    end
end

% Map the smoothed areas to the ATLAS shape
color = zeros(size(faceA, 1), 1);
color(IdxCUT) = newarea(IdxMatch);

%% Visualization of Results

figure;

% Plot live-shape with smoothed area values
subplot(1, 2, 1);
plotmesh(node, face, 'FaceVertexCData', newarea, 'EdgeColor', 'none');
caxis([min(newarea) max(newarea)]);
view([180 0]);
title('Live-Shape');
axis off;

% Plot ATLAS with corresponding area values
subplot(1, 2, 2);
plotmesh(nodeA, faceA, 'FaceVertexCData', color, 'EdgeColor', 'none');
caxis([min(newarea) max(newarea)]);
view([180 0]);
title('ATLAS');
axis off;
