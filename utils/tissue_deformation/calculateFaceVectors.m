% Function to calculate face vectors
function [tv, rv, sv] = calculateFaceVectors(faceNodes)
    % Compute vector from first to second node
    rv = faceNodes(2, :) - faceNodes(1, :);
    rv = rv / norm(rv);  % Normalize rv

    % Compute the normal vector (tv) using cross product
    tv = cross(rv, faceNodes(3, :) - faceNodes(1, :));
    tv = tv / norm(tv);  % Normalize tv

    % Compute the perpendicular vector (sv) using cross product
    sv = cross(rv, tv);
    sv = sv / norm(sv);  % Normalize sv
end