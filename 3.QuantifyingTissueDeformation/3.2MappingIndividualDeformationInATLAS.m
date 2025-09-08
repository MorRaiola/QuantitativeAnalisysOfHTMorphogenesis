%% 3.2 Plotting Deformation onto the ATLAS Shape using Face-to-Face Matching
clear; close all; clc;

% Read data from the Excel file
Excel = xlsread('\Code_data\2. IntegratingMultipleLiveImagesIntoAConsensusTemporalReference\EmbryoStage.xlsx'); 

% Loop through each row in the Excel file
for cla = 1:size(Excel, 1)
    
    % Define base folder and relevant paths
    baseFolder = '\Code_data\2. IntegratingMultipleLiveImagesIntoAConsensusTemporalReference';
    clusterID = num2str(Excel(cla, 1));
    embryoID = num2str(Excel(cla, 2));
    
    deformationFolder = fullfile(baseFolder, ['Cluster', clusterID], 'Deformation');
    atlasData = fullfile(baseFolder, 'Atlas', ['remesh', clusterID, '.ply']);
    outputFolder = fullfile('C:\Users\mraiola\Desktop\im\', ['Gr', clusterID], 'ATLAS\');
    
    % Create output directory if it doesn't exist
    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end
    
    % Load and reorder PLY file data
    [vertices, faces] = read_ply(atlasData);
    vertices = vertices(:, [3, 1, 2]);  % Reorder vertices
    faces = faces(:, [3, 1, 2]);        % Reorder face indices
    
    myo.vertices = vertices;
    myo.faces = faces;
    
    %% MAP GROWTH RATE
    % Load deformation and corresponding face indices
    load(fullfile(deformationFolder, ['JCCsmooth', embryoID, '.mat']));
    load(fullfile(baseFolder, 'Mapping', ['Cluster', clusterID], 'VertexCorrespondence', ['Embryo', embryoID], 'IdxCUT.mat'));
    load(fullfile(baseFolder, 'Mapping', ['Cluster', clusterID], 'VertexCorrespondence', ['Embryo', embryoID], 'IdxMatch.mat'));
    
    % Map the growth rate to the mesh
    JCC = newJCC;
    color = zeros(size(faces, 1), size(JCC, 2));
    for i = 1:size(color, 2)
        color(IdxCUT, i) = JCC(IdxMatch, i);
    end
    
    % Convert data to RGB and write the mesh to a PLY file
    RGB = dataToRGB(color, 1);
    outputFile = fullfile(outputFolder, [embryoID, '_GrowthRate.ply']);
    delete(outputFile);
    writeMesh_plyModify(outputFile, myo.vertices, myo.faces, color, RGB);

    %% MAP ANISOTROPY
    % Load anisotropy data and corresponding face indices
    load(fullfile(deformationFolder, ['Anisotropysmooth', embryoID, '.mat']));
    load(fullfile(baseFolder, 'Mapping', ['Cluster', clusterID], 'VertexCorrespondence', ['Embryo', embryoID], 'IdxCUT.mat'));
    load(fullfile(baseFolder, 'Mapping', ['Cluster', clusterID], 'VertexCorrespondence', ['Embryo', embryoID], 'IdxMatch.mat'));
    
    % Map the anisotropy to the mesh
    Anisotropy = newAnisotropy;
    color = zeros(size(faces, 1), size(Anisotropy, 2));
    for i = 1:size(color, 2)
        color(IdxCUT, i) = Anisotropy(IdxMatch, i);
    end
    
    % Convert data to RGB and write the mesh to a PLY file
    RGB = dataToRGB(color, 1);
    outputFile = fullfile(outputFolder, [embryoID, '_Anisotropy.ply']);
    delete(outputFile);
    writeMesh_plyModify(outputFile, myo.vertices, myo.faces, color, RGB);
    
    % Clear variables except for Excel
    clearvars -except Excel;
end

%% Defining mean and std for growth rate and anysotropy ∀ Gr
clear; close all; clc;

% Define cluster identifiers
clusters = [2, 3, 4, 5, 6, 7, 8, 9];

% Base folder path for results
baseFolder = '\Code_data\2. IntegratingMultipleLiveImagesIntoAConsensusTemporalReferencea';
outputBaseFolder = 'V:\CellReportsMethods_2024\Data\Figure4_Quantifying_Tissue_Deformation_During_Early_Morphogenesis';

% Loop through each cluster
for clustIdx = 1:numel(clusters)
    clusterID = clusters(clustIdx);
    
    % Define file paths and names
    dataFile = fullfile('\Source_Data\Figure 3\3A', sprintf('Gr%d.ply', clusterID));
    deformationFolder = fullfile(baseFolder, sprintf('Cluster%d', clusterID), 'Deformation');
    outputFolder = fullfile(outputBaseFolder, sprintf('Gr%d', clusterID), 'StepWise');
   
    % Read PLY file data
    [vertices, faces] = read_ply(dataFile);
    myo.vertices = vertices;
    myo.faces = faces;
    
    % Process both JCCsmooth and Anisotropysmooth files
    for typeIdx = 1:2
        if typeIdx == 1
            fileType = 'JCCsmooth';
            resultPrefix = 'GrowthRate';
            resultData = 'newJCC';
        else
            fileType = 'Anisotropysmooth';
            resultPrefix = 'Anisotropy';
            resultData = 'newAnisotropy';
        end
        
        % Get list of files
        filePattern = sprintf('*%s*', fileType);
        files = dir(fullfile(deformationFolder, filePattern));
        fileNames = {files.name};
        colorData = zeros(size(faces, 1), numel(fileNames));
        
        % Process each file
        for fileIdx = 1:numel(fileNames)
            filePath = fullfile(deformationFolder, fileNames{fileIdx});
            load(filePath);
            
            % Extract embryo number from file name
            embryoNum = regexp(files(fileIdx).name, '\d+', 'match');
            embryoNum = embryoNum{1};
            
            % Load vertex correspondence files
            idxCutFile = fullfile(baseFolder, sprintf('Cluster%d', clusterID), 'VertexCorrespondence', sprintf('Embryo%s', embryoNum), 'IdxCUT.mat');
            idxMatchFile = fullfile(baseFolder, sprintf('Cluster%d', clusterID), 'VertexCorrespondence', sprintf('Embryo%s', embryoNum), 'IdxMatch.mat');
            load(idxCutFile);
            load(idxMatchFile);
            
            % Update color data based on current result type
            values = eval(resultData);
            colorData(IdxCUT, fileIdx) = values(IdxMatch, 1);
        end
        
        % Handle zero values and compute statistics
        colorData(any(colorData == 0, 2), :) = 0;
        if size(colorData, 2) > 1
            colorMean = mean(colorData, 2);
            colorSTD = std(colorData, [], 2);
        else
            colorMean = colorData;
            colorSTD = zeros(size(colorMean, 1), 1);
        end
        
        % Prepare and create output folder if necessary
        if ~exist(outputFolder, 'dir')
            mkdir(outputFolder);
        end
        
        % Save mean and standard deviation data
        for statIdx = 1:2
            if statIdx == 1
                data = colorMean;
                suffix = 'Mean';
            else
                data = colorSTD;
                suffix = 'Std';
            end
            
            if sum(data)>1
                rgbData = dataToRGB(data,  1, typeIdx,suffix);
                filePath = fullfile(outputFolder, sprintf('%s_%s.ply', suffix, resultPrefix));
                if exist(filePath, 'file')
                    delete(filePath);
                end
                writeMesh_plyModify(filePath, myo.vertices, myo.faces, data, rgbData);
            end
        end

    end
    
    % Clear variables for the next iteration
    clearvars -except clusters baseFolder outputBaseFolder
end
