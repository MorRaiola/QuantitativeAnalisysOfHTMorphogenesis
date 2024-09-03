%% 3.1 Individual Tissue Deformation (deformation betweeen consecutive stages)
% Live-shape deformation (JCC,Anisotropy)
% Individual Tissue Deformation Analysis
% This script calculates tissue deformation between consecutive stages.

clear; close all; clc;

% Load the Excel file containing embryo stage information
excelFilePath = '\\tierra.cnic.es\SC\LAB_MT\RESULTADOS\Morena\EmbryoStage.xlsx';
excelData = xlsread(excelFilePath);

% Process each entry in the Excel file
for i = 1:size(excelData, 1)
    clusterId = excelData(i, 1);
    embryoId = excelData(i, 2);
    timeStart = excelData(i, 3);
    timeEnd = excelData(i, 4);

    % Define file paths
    deformationFolder = (['\\tierra.cnic.es\SC\LAB_MT\RESULTADOS\Morena\Mapping\Cluster' num2str(clusterId) '\Deformation']);
    dataFolder = fullfile(['\\tierra.cnic.es\SC\LAB_MT\RESULTADOS\Morena\Embryos\Embryo' num2str(embryoId) '\Shapes\CC']);

    % Load the initial and final node/face data
    [node1, face1] = loadMeshData(dataFolder, timeStart);
    [node2, face2] = loadMeshData(dataFolder, timeEnd);

    % Calculate the local coordinate systems and deformation metrics
    [F, JCC, anisotropy, globalAutovet, barycenter2] = calculateDeformationMetrics(node1, face1, node2, face2);

    % Save the results
    saveDeformationData(deformationFolder, embryoId, JCC, anisotropy, globalAutovet, barycenter2);

    % Clear variables for the next iteration
    clearvars -except excelData
end


