%% 3.1 Individual Tissue Deformation (deformation betweeen consecutive stages)
% Live-shape deformation (JCC,Anisotropy)
% Individual Tissue Deformation Analysis
% This script calculates tissue deformation between consecutive stages.

clear; close all; clc;

% Load the Excel file containing embryo stage information
excelFilePath = '\Code_data\2. IntegratingMultipleLiveImagesIntoAConsensusTemporalReference\EmbryoStage.xlsx';
excelData = xlsread(excelFilePath);
clusterFolder = '\Code_data\2. IntegratingMultipleLiveImagesIntoAConsensusTemporalReference\Cluster';
embryoFolder = '\Code_data\1. EstimatingIndividualLiveImageMotion\Embryo';

% Process each entry in the Excel file
for i = 1:size(excelData, 1)
    clusterId = excelData(i, 1);
    embryoId = excelData(i, 2);
    timeStart = excelData(i, 3);
    timeEnd = excelData(i, 4);

    % Define file paths
    deformationFolder = ([clusterFolder num2str(clusterId) '\Deformation']);
    dataFolder = fullfile([embryoFolder num2str(embryoId) '\Shapes\CC']);

    % Load the initial and final node/face data
    [node1, face1] = loadMeshData(dataFolder, timeStart);
    [node2, face2] = loadMeshData(dataFolder, timeEnd);

    % Calculate the local coordinate systems and deformation metrics
    [F, JCC, anisotropy, globalAutovet, barycenter2] = calculateDeformationMetrics(node1, face1, node2, face2);
    figure
    plotmesh(node2,face2,'FaceVertexCData',anisotropy,'EdgeColor','none')
    % Save the results
    saveDeformationData(deformationFolder, embryoId, JCC, anisotropy, globalAutovet, barycenter2);

    % Clear variables for the next iteration
    clearvars -except excelData
end


%% Smooth deformation in 20-pixel distance
clear; close all; clc;

Excel = xlsread('\Code_data\2. IntegratingMultipleLiveImagesIntoAConsensusTemporalReference\EmbryoStage.xlsx');
clusterFolder = '\Code_data\2. IntegratingMultipleLiveImagesIntoAConsensusTemporalReference\Cluster';

for cla = 1:size(Excel, 1)
    % File paths and variables
    clusterPath = [clusterFolder num2str(Excel(cla, 1))];
    DataCC = fullfile(clusterPath, '3D', 'Embryo', [num2str(Excel(cla, 2)) '.ply']);
    Folder = fullfile(clusterPath, 'Deformation');
    embryo = Excel(cla, 2);
    
    % Load and preprocess mesh
    [node, face] = read_ply(DataCC);
    node = node(:, [3 1 2]);
    face = face(:, [3 1 2]);
    centroid = meshcentroid(node, face);
    
    % Set up geodesic library
    global geodesic_library;                
    geodesic_library = 'geodesic_release';
    
    % Initialize geodesic mesh and algorithm
    tic;
    mesh = geodesic_new_mesh(node, face);
    algorithm = geodesic_new_algorithm(mesh, 'exact');
    
    % Compute indices
    index = cell(size(face, 1), 1);
    for i = 1:size(face, 1)
        vertex_id = knnsearch(node, centroid(i, :));
        source_points = {geodesic_create_surface_point('vertex', vertex_id, node(vertex_id, :))};
        geodesic_propagate(algorithm, source_points, [], 20);
        
        [~, distances] = geodesic_distance_and_source(algorithm);
        dist = find(distances < 1e10);
        
        index{i} = all(ismember(face, dist), 2);
    end
    toc;
    
    % Save result
    save(fullfile(Folder, ['indexSmooth' num2str(embryo) '.mat']), 'index');
end


%% Save in live-shape smoothed deformation map
clear; close all; clc;

% Load the Excel data
excelFile = '\Code_data\2. IntegratingMultipleLiveImagesIntoAConsensusTemporalReference\EmbryoStage.xlsx';
excelData = xlsread(excelFile); 

% Define base folder paths
baseFolder = '\Code_data\2. IntegratingMultipleLiveImagesIntoAConsensusTemporalReference';
outputBaseFolder = 'C:\Users\mraiola\Desktop\im\Gr';

% Loop through each row in the Excel file
for i = 1:size(excelData, 1)
    % Extract Gr and embryo information
    clusterNum = num2str(excelData(i, 1));
    embryoNum = num2str(excelData(i, 2));
    
    % Define folders and file paths
    clusterFolder = fullfile(baseFolder, ['Cluster' clusterNum], 'Deformation');
    dataFile = fullfile(baseFolder, ['Cluster' clusterNum], '3D', 'Embryo', [embryoNum '.ply']);
    outputFolder = fullfile(outputBaseFolder, clusterNum, 'Individual');
    
    % Create output folder if it doesn't exist
    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end
    
    % Load mesh data
    [vertices, faces] = read_ply(dataFile);
    vertices = vertices(:, [3 1 2]);
    faces = faces(:, [3 1 2]);

    % Load deformation data
    baricentroFile = fullfile(clusterFolder, ['baricentro2' embryoNum '.mat']);
    autovetFile = fullfile(clusterFolder, ['GlobalautovetCC' embryoNum '.mat']);
    growthRateFile = fullfile(clusterFolder, ['JCCsmooth' embryoNum '.mat']);
    
    load(baricentroFile, 'baricentro2');
    load(autovetFile, 'autovettoriGlob');
    load(growthRateFile, 'newJCC');

    % Generate color map and write GrowthRate mesh
    rgbGrowthRate = dataToRGB(newJCC, 1);
    growthRateOutput = fullfile(outputFolder, [embryoNum '_GrowthRate.ply']);
    
    if exist(growthRateOutput, 'file')
        delete(growthRateOutput);
    end
    
    writeMesh_plyModify(growthRateOutput, vertices, faces, newJCC, rgbGrowthRate);
    
    % Load and process Anisotropy data
    anisotropyFile = fullfile(clusterFolder, ['Anisotropysmooth' embryoNum '.mat']);
    load(anisotropyFile, 'newAnisotropy');
    
    rgbAnisotropy = dataToRGB(newAnisotropy, 1);
    anisotropyOutput = fullfile(outputFolder, [embryoNum '_Anisotropy.ply']);
    
    if exist(anisotropyOutput, 'file')
        delete(anisotropyOutput);
    end
    
    writeMesh_plyModify(anisotropyOutput, vertices, faces, newAnisotropy, rgbAnisotropy);
    
    % Clear variables except for the main Excel data
    clearvars -except excelData baseFolder outputBaseFolder;
end



