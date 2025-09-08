%% 1.4 Continuous Description of Heart Tube (HT) Morphogenesis (Mesh)
clear; close all; clc;

% Add MIJ and ImageJ Java paths
addpath('C:\Users\mraiola\Downloads\fiji-win64\Fiji.app\scripts');
javaaddpath('C:\Program Files\MATLAB\R2020a\java\jar\mij.jar'); % MIJ Java folder
javaaddpath('C:\Program Files\MATLAB\R2020a\java\jar\ij-1.52i.jar'); % ImageJ JAR file

% Load the time-lapse image data
Embryo = '\Code_Data\1. EstimatingIndividualLiveImageMotion\Embryo1\Embryo1\Data\Embryo1.tif';
transformationPath = '\Code_Data\1. EstimatingIndividualLiveImageMotion\Embryo1\disp.mat';
outputMeshPath = '\Code_Data\1. EstimatingIndividualLiveImageMotion\Embryo1\Shapes\CC';

% Load the image data
imp = ij.IJ.openImage(Embryo); % Load image (ensure RGB images are converted to grayscale for feature extraction)
Im = ImagePlus2array(imp);
Im = squeeze(Im);

% Load transformation data
load(transformationPath, 'disp'); % Load displacement data
timepoints = size(disp, 1) + 1;
midline = uint8(timepoints / 2);

% Generate 3D mesh from the segmented image
inputTiffPath = fullfile(baseFolder, ['Embryo' num2str(embryoID)], 'Segmentation', ['Embryo' num2str(embryoID) '.tif']);
outputPlyPath = fullfile(baseFolder, ['Embryo' num2str(embryoID)], 'Shapes', ['CC' num2str(midline) '.ply']);
[newVertices, newFaces] = Mesh3DFromImage(inputTiffPath, outputPlyPath, false,false);

% Initialize surface structure for the heart tube
srfCC.vertices = newVertices;
srfCC.faces = newFaces;
srfCC.disps = zeros(size(srfCC.vertices));

% Initialize original surface structure
srfCC0 = srfCC;

% Set origin and spacing for 3D morphing
orig = struct('x', 0, 'y', 0, 'z', 0);
spacing = struct('x', 1, 'y', 1, 'z', 1);

% Create the motion profile for the heart tube

tic
% Forward morphing: Midline to first time point
for t = 1:midline-1
    newsrf = morphSurf3D(srfCC, disp{midline-t, 1}, orig, spacing, Im);
    srfCC.vertices = newsrf.vertices;
    write_ply(srfCC.vertices, srfCC.faces, [outputMeshPath num2str(midline-t) '.ply']);
end

% Reset surface structure for backward morphing
srfCC.vertices = newVertices;
srfCC.faces = newFaces;

% Backward morphing: Midline to last time point
for t = midline:timepoints-1
    newsrf = morphSurf3D(srfCC, disp{t, 1}, orig, spacing, Im);
    srfCC.vertices = newsrf.vertices;
    write_ply(srfCC.vertices, srfCC.faces, [outputMeshPath num2str(midline) '.ply']);
end
toc

