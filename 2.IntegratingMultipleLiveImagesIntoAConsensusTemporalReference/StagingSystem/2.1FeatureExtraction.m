%% Morphometric Feature Definition
% This script computes a morphometric feature ratio (height/width) for a collection
% of 3D meshes stored in .ply files. The script loads each mesh, allows the user to 
% manually select points corresponding to width and height, and then calculates the 
% ratio of these distances for each mesh.

clear all; close all; clc;

% Directory containing the .ply files
InDataDir = '\Source_Data\Figure 3\3A';

% Get a list of all .ply files in the directory
inFiles = dir(fullfile(InDataDir, '*.ply'));
inNames = cell(1, length(inFiles));

% Store the full paths of the .ply files in inNames
for i = 1:length(inFiles)
    inNames{i} = fullfile(InDataDir, inFiles(i).name);
end

% Initialize arrays to store the width, height, and computed feature values
w = zeros(length(inFiles), 2);
h = zeros(length(inFiles), 2);
feature = zeros(length(inFiles), 1);

% Loop through each .ply file to compute the morphometric feature
for i = 1:length(inNames)
    % Load the 3D mesh data
    [node, face] = read_ply(inNames{i});
    
    % Plot the mesh for manual point selection
    plotmesh(node, face);
    
    % Manually select points on the mesh corresponding to width and height
    disp('Select two points for width (w1, w2), then press Enter.');
    w(i,:) = [cursor_info(3).Position cursor_info(4).Position];
    h(i,:) = [cursor_info(1).Position cursor_info(2).Position];
 
    % Calculate the feature ratio (height/width) for the current mesh
    feature(i) = norm(h(i,1:3) - h(i,4:6)) / norm(w(i,1:3) - w(i,4:6));
end

% Display the computed features
disp('Computed morphometric feature ratios (height/width):');
disp(feature);
