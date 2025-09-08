%% 1.5 Validating the Motion Estimation
clear; close all; clc;

% Add MIJ Java paths for MATLAB to interact with ImageJ
addpath('C:\Users\mraiola\Downloads\fiji-win64\Fiji.app\scripts');
javaaddpath('C:\Program Files\MATLAB\R2020a\java\jar\mij.jar');
javaaddpath('C:\Program Files\MATLAB\R2020a\java\jar\ij-1.52i.jar');

% Define the base folder and embryo ID
baseFolder = '\Source_Data\Figure 2\2F\25%\Validating the Motion Estimation\';
imageFolder = '\Code_Data\1. EstimatingIndividualLiveImageMotion\';
embryoID = 1;

% Load transformation data
load(fullfile(baseFolder, 'disp.mat'), 'disp');

% Load manual tracking data from Imaris
trackingFile = fullfile(baseFolder, 'Tracking', ['e' num2str(embryoID) '_Statistics'], 'GroundTruth.csv');
Excel = xlsread(trackingFile);

% Define parameters
downscaleFactor = 0.25; % Downscaling factor
resImage = 0.593; % Isotropic voxel size in micrometers (µm)
voxelSize = resImage / downscaleFactor; % Factor to convert to micrometers
cellDiameter = 20; % Diameter of a cardiomyocyte in µm

% Load the image for validation
imageFile = fullfile(imageFolder, ['Embryo' num2str(embryoID)], 'Data, ['Embryo' num2str(embryoID) '.tif']);
imp = ij.IJ.openImage(imageFile); 
Im = squeeze(ImagePlus2array(imp));
Im = Im(:,:,:,1); % Use the first time point

% Initialize matrices for tracking data
timepoints = size(unique(Excel(:,7)), 1);
ntracks = length(find(Excel(:,7) == timepoints));
d = 1e9; % Large constant to ensure uniqueness
index = unique(mod(Excel(:,8), d)); 
trackId = zeros(timepoints, ntracks);
for i = 1:length(index)
    trackId(:,i) = find(mod(Excel(:,8), d) == index(i));
end

% Set origin and spacing for morphing operations
orig = struct('x', 0, 'y', 0, 'z', 0);
spacing = struct('x', 1, 'y', 1, 'z', 1);

% Initialize structures for validation
prop = Excel;
punc = Excel;
midline = uint8(timepoints / 2);
srfCC.vertices = Excel(trackId(midline,:), 1:3);

% Sequential error test set creation
tic
for t = 1:midline-1
    newsrf = morphSurf3D(srfCC, disp{midline-t, 1}, orig, spacing, Im(:,:,:));
    srfCC.vertices = newsrf.vertices;
    prop(trackId(midline-t,:), 1:3) = newsrf.vertices;
end   

srfCC.vertices = Excel(trackId(midline,:), 1:3);
for t = midline:timepoints-1
    newsrf = morphSurf3D(srfCC, disp{t, 1}, orig, spacing, Im(:,:,:));
    srfCC.vertices = newsrf.vertices;
    prop(trackId(t+1,:), 1:3) = newsrf.vertices;
end   
toc

% Save the results to CSV
outputTestFile = fullfile(baseFolder, 'Tracking', ['e' num2str(embryoID) '_Statistics'], 'Testset.csv');
writetable(prop, outputTestFile);

% Punctual error test set creation
tic
for t = 1:midline-1
    srfCC.vertices = Excel(trackId(midline-t+1,:), 1:3);
    newsrf = morphSurf3D(srfCC, disp{midline-t, 1}, orig, spacing, Im(:,:,:));
    punc(trackId(midline-t,:), 1:3) = newsrf.vertices;
end   

for t = midline:timepoints-1
    srfCC.vertices = Excel(trackId(midline,:), 1:3);
    newsrf = morphSurf3D(srfCC, disp{t, 1}, orig, spacing, Im(:,:,:));
    punc(trackId(t+1,:), 1:3) = newsrf.vertices;
end   
toc

% Evaluate the error using Eulerian distance
GroundTruth = Excel;
Testset = prop; % Alternatively, use punc for punctual error

EulDist = zeros(timepoints, ntracks);
for i = 1:timepoints
    for t = 1:size(trackId, 2)
        EulDist(i, t) = sqrt(sum((GroundTruth(trackId(i,t), 1:3) - Testset(trackId(i,t), 1:3)).^2));
    end
end

EulDist = EulDist * voxelSize; % Convert error to micrometers
EulDist(midline,:) = []; % Exclude midline data
meanError = mean(EulDist(:), 'all', 'omitnan');
maxError = max(EulDist(:), [], 'omitnan');
stdError = std(EulDist(:), 'omitnan');

% Percentage of cells with errors within specific thresholds
below10um = size(find(EulDist <= (cellDiameter / 2)), 1) * 100 / (size(EulDist, 2) * size(EulDist, 1));
between10umAnd20um = size(find(EulDist > (cellDiameter / 2) & EulDist <= cellDiameter), 1) * 100 / (size(EulDist, 2) * size(EulDist, 1));
between20umAnd40um = size(find(EulDist > cellDiameter & EulDist <= (2 * cellDiameter)), 1) * 100 / (size(EulDist, 2) * size(EulDist, 1));
above40um = size(find(EulDist > (2 * cellDiameter)), 1) * 100 / (size(EulDist, 2) * size(EulDist, 1));

% Identify and isolate outliers
[TF, L, U, Center] = isoutlier(EulDist(:)); % Find outliers
[x, y] = find(EulDist > U);
outlierData = Testset(trackId(:, unique(y)), :);

