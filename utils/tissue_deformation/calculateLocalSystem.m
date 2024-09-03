% Function to calculate the local coordinate system and barycenters
function [localSystem, barycenter] = calculateLocalSystem(nodes, faces)
    numFaces = size(faces, 1);
    localSystem = zeros(numFaces, 9);
    barycenter = zeros(numFaces, 3);

    for z = 1:numFaces
        barycenter(z, :) = mean(nodes(faces(z, 1:3), :), 1);
        [tv, rv, sv] = calculateFaceVectors(nodes(faces(z, 1:3), :));
        localSystem(z, :) = [rv, sv, tv];
    end
end
