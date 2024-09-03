function [F, JCC, anisotropy, globalAutovet, barycenter2] = calculateDeformationMetrics(node1, face1, node2, face2)
    % Calculate local coordinate systems and barycenters for both sets of nodes and faces
    [localSystem1, barycenter1] = calculateLocalSystem(node1, face1);
    [localSystem2, barycenter2] = calculateLocalSystem(node2, face2);

    numFaces = size(face1, 1); % Number of faces
    
    % Initialize output variables
    F = cell(numFaces, 1);
    JCC = zeros(numFaces, 1);
    strainMax = zeros(numFaces, 1);
    anisotropy = zeros(numFaces, 1);
    globalAutovet = cell(numFaces, 1);

    % Fixed vectors for transformation
    yv = [0 1 0]; % Longitudinal vector
    xv = [1 0 0]; % Circumferential vector
    zv = [0 0 1]; % Normal vector

    % Loop through each face to compute deformation metrics
    for z = 1:numFaces
        % Extract nodes corresponding to the current face
        faceNodes1 = node1(face1(z, 1:3), :);
        faceNodes2 = node2(face2(z, 1:3), :);

        % Compute deformation metrics for the current face
        [F{z}, JCC(z), anisotropy(z), globalAutovet{z}] = ...
            calculateFaceDeformation(faceNodes1, faceNodes2, localSystem1(z, :), localSystem2(z, :), xv, yv, zv);
    end
end
