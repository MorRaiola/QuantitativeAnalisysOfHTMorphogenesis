%% 1.6 Validating HT Morphogenesis Description
% From segmentation to mesh of the GroundTruth
clear; close all; clc;

% Add MIJ and ImageJ Java paths
addpath('C:\Users\mraiola\Downloads\fiji-win64\Fiji.app\scripts');
javaaddpath('C:\Program Files\MATLAB\R2020a\java\jar\mij.jar');  % MIJ Java folder
javaaddpath('C:\Program Files\MATLAB\R2020a\java\jar\ij-1.52i.jar');  % ImageJ JAR file

% Load cell segmentation data
embryoFile = '\1. EstimatingIndividualLiveImageMotion\1.6Error2\Cells\Reslice of thalf.tif';
output = '\1. EstimatingIndividualLiveImageMotion\1.6Error2\Cells\GT'
imp = ij.IJ.openImage(embryoFile);
Seg = squeeze(ImagePlus2array(imp));
Im = Seg(:,:,:,1);  % Use the first time point

% Create mesh for each segmented cell
uniqueCells = unique(Im);
for i = 1:length(uniqueCells) - 1
    IM = Im;
    IM(IM ~= uniqueCells(i)) = 0;
    IM(IM == uniqueCells(i)) = 1;  % Select single cell by label value
    cleancc = fillholes3d(logical(IM > 0), 1);
    [node, elem, face] = v2m(cleancc, 0.1, 0.5, 0.5);  % Generate mesh

    % Smooth the mesh
    conn = meshconn(face(:, 1:3), size(node, 1));
    n1{1} = node;
    niter = 6;
    for iter = 2:niter + 1
        n1{iter} = smoothsurf(n1{iter-1}, [], conn, 1, 0.9, 'laplacian');
    end
    node = n1{end}(:, [2, 1, 3]);
    [newnode2, newface2] = surfreorient(node, face);

    % Save the mesh data
    figure;
    plotmesh(newnode2, newface2, 'facealpha', 0.7);
    saveFile = fullfile(output, ['tmezzi0000' num2str(i) '.ply']);
    write_ply(newnode2, newface2, saveFile);
end

%% Continuous Description of Cell Shapes during Morphogenesis
clear; close all; clc;

% Load final time point segmentation
Im = loadtiff('\1. EstimatingIndividualLiveImageMotion\1.6Error2\Cells\Reslice of tend.tif');
load('\1. EstimatingIndividualLiveImageMotion\Embryo2\disp.mat');
output = '\1. EstimatingIndividualLiveImageMotion\1.6Error2\Cells\Test';

% Morph each cell's mesh over time
for i = 1:length(unique(Im)) - 1
    [newnode2, newface2] = read_ply([output '\tmezzi0000' num2str(i) '.ply']);

    srfCC.vertices = newnode2;
    srfCC.faces = newface2;
    srfCC.disps = zeros(size(srfCC.vertices));  % Initialize displacement

    orig.x = 0; orig.y = 0; orig.z = 0;
    spacing.x = 1; spacing.y = 1; spacing.z = 1;

    timepoints = size(disp, 1) + 1;
    midline = uint8(timepoints / 2);

    % Morph mesh backward in time
    for t = 1:midline - 1
        newsrf = morphSurf3D(srfCC, disp{midline - t, 1}, orig, spacing, Im(:,:,:));
        srfCC.vertices = newsrf.vertices;
    end
    saveFilePrimo = fullfile(output, ['primo0000' num2str(i) '.ply']);
    write_ply(srfCC.vertices, srfCC.faces, saveFilePrimo);

    srfCC.vertices = newnode2;  % Reset to original vertices

    % Morph mesh forward in time
    for t = midline:timepoints - 1
        newsrf = morphSurf3D(srfCC, disp{t, 1}, orig, spacing, Im(:,:,:));
        srfCC.vertices = newsrf.vertices;
    end
    saveFileUltimo = fullfile(output', ['ultimo0000' num2str(i) '.ply']);
    write_ply(srfCC.vertices, srfCC.faces, saveFileUltimo);

    clearvars -except Im disp i;  % Clear unnecessary variables for the next iteration
end

%% Error Calculation: Cosine Similarity between Tissue Deformation and Cell Division Direction
clear; close all; clc;

% Load initial mesh and segmentation
cellMid = loadtiff('\1. EstimatingIndividualLiveImageMotion\1.6Error2\Cells\Reslice of t1.tif');
[V, F] = read_ply('\1. EstimatingIndividualLiveImageMotion\Embryo2\Shape\CC1.ply');
Cs = zeros(length(unique(Im)) - 1, 1); 
angles = zeros(length(unique(Im)) - 1, 1); 

% Display initial mesh
plotmesh(V, F, 'facealpha', 0.2, 'edgecolor', 'none');
hold on;

% Calculate cosine similarity for each cell
for i = 1:length(unique(Im)) - 1
    % Load cell mesh from the test set (_primo)
    [V_test, F_test] = read_ply(['\1. EstimatingIndividualLiveImageMotion\1.6Error2\Cells\Test\primo0000' num2str(i) '.ply']);
    line_test = fitLine3d(V_test);
    
    % Load corresponding cell mesh from the Ground Truth (_primo)
    [V_gt, F_gt] = read_ply(['1. EstimatingIndividualLiveImageMotion\1.6Error2\Cells\GT\primo' num2str(i) '.ply']);
    line_gt = fitLine3d(V_gt);

    % Calculate cosine similarity between test set and Ground Truth for _primo
    Cs(i, 1) = getCosineSimilarity(line_test(1, 4:6), line_gt(1, 4:6));

    % Load cell mesh from the test set (_ultimo)
    [V_test, F_test] = read_ply(['\1. EstimatingIndividualLiveImageMotion\1.6Error2\Cells\Test\ultimo0000' num2str(i) '.ply']);
    line_test = fitLine3d(V_test);
    
    % Load corresponding cell mesh from the Ground Truth (_ultimo)
    [V_gt, F_gt] = read_ply(['1. EstimatingIndividualLiveImageMotion\1.6Error2\Cells\GT\ultimo' num2str(i) '.ply']);
    line_gt = fitLine3d(V_gt);

    % Calculate cosine similarity between test set and Ground Truth for _ultimo
    Cs(i, 2) = getCosineSimilarity(line_test(1, 4:6), line_gt(1, 4:6));
end

