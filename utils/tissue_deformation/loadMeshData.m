% Function to load mesh data
function [nodes, faces] = loadMeshData(folderPath, timePoint)
    filePath = ([folderPath  num2str(timePoint) '.ply']);
    [nodes, faces] = read_ply(filePath);
    nodes = nodes(:, [3, 1, 2]);  % Reorder to match required format
    faces = faces(:, [3, 1, 2]);
end
