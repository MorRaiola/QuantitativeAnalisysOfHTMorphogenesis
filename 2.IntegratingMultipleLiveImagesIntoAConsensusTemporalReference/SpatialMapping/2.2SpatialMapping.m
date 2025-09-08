%% Mapping live-shapes into the relative ATLAS Gr shape
%% 2.2.1 Rigid registration atlas-->live images with TGMM
clear; close all; clc;
Excel = xlsread('\\2. IntegratingMultipleLiveImagesIntoAConsensusTemporalReference\EmbryoStage.xlsx'); 

for cla = 1:size(Excel, 1)
    % Load file paths
    AtlasFolder = sprintf('\Source_Data\Figure 3\3A\Gr%d.ply', Excel(cla, 1));
    embryo = Excel(cla, 2);
    EmbryoFolder = sprintf('\\2. IntegratingMultipleLiveImagesIntoAConsensusTemporalReference\\Cluster%d\\3D\\Embryo\\%d.ply', Excel(cla, 1), embryo);
    OutFolder = sprintf('\\2. IntegratingMultipleLiveImagesIntoAConsensusTemporalReferencen\\Cluster%d\\3D\\Atlas\\%d.ply', Excel(cla, 1), embryo);

    [nodeA, faceA] = read_ply(AtlasFolder); % Read PLY file ATLAS
    [nodeE, faceE] = read_ply(EmbryoFolder); % Read PLY file Embryo 
    Tset{1,1}=nodeA;
    % Tset{1,2}=nodeE;
    % Apply TGMM https://es.mathworks.com/matlabcentral/fileexchange/63693-robust-group-wise-registration-of-point-sets-using-multi-resolution-t-mixture-model
    [MU, Transform, TrainingSet, UP, PP, Mcoeffs, nu, convg, SSM] = TMMgroupwiseReg(Tset, uint16(size(nodeE, 1)/ 2), 150, 1, nodeE);

    write_ply(TrainingSet.TransfPts, faceA, OutFolder);
    save(sprintf('\\2. IntegratingMultipleLiveImagesIntoAConsensusTemporalReference\\Cluster%d\\Disp\\Rigid\\%d.mat', Excel(cla, 1), embryo), 'Transform');
end

% Rotate again the data in the ATLAS reference
tempSR = bsxfun(@times, Transform.R, Transform.s);
trPts = tempSR * TrainingSet.TransfPts';
Xo = bsxfun(@plus, trPts', Transform.t);

%% 2.2.2 Manual Cut of missing parts in MeshLab --> Cut_(embryo).ply / Fillgaps in Meshlab --> Cut_(embryo).1.ply

All the manual Cut are in '\Supplementary_Data\Figure 3\3E\SpatialMapping\Grx\Atlas'

%% 2.2.3 Mask ATLAS_Cut, creating the ATLAS mask without the missing parts (IFTs or OFT)
clear; close all; clc;

Excel = xlsread('\\2. IntegratingMultipleLiveImagesIntoAConsensusTemporalReference\EmbryoStage.xlsx'); 

for cla = 1:size(Excel, 1)
    Folder = sprintf('\\2. IntegratingMultipleLiveImagesIntoAConsensusTemporalReference\\Cluster%d', Excel(cla, 1));
    embryo = Excel(cla, 2);
    time = Excel(cla, 4);
    OutFolder = sprintf('%s\\Images\\Atlas\\%d.tif', Folder, embryo);

    [node, face] = read_ply(sprintf('%s\\3D\\Atlas\\Cut_%d.1.ply', Folder, embryo)); % Filled Cut mesh
    SegTest = loadtiff(sprintf('%s\\Images\\Embryo\\%d.%d.tif', Folder, embryo, time)); % Embryo segmentation as size reference

    % Create the segmentation from surface mesh
    [dx, dy, dz] = deal(1:size(SegTest, 1), 1:size(SegTest, 2), 1:size(SegTest, 3));
    imag = surf2volz(node, face, dx, dy, dz); % From mesh to binary mask
    options.overwrite = true;
    saveastiff(imag, OutFolder, options);

    clearvars -except Excel;
end



%% 2.2.4 Import the segmentation in Fiji. Fill gap manually + close holes --> Fill_(embryo).tif



%% 2.2.4 Import the segmentation in Fiji. Fill gap manually + close holes --> Fill_(embryo).tif5 Mask live-shape
clear; close all; clc;

Folder = '\\Supplementary_Data\Figure 2\2I\Segmenting Heart Tissue';
Folderdisp ='\\1. EstimatingIndividualLiveImageMotion';

embryo =1;
Embryo = sprintf('%s\\e0%d.tif', Folder, embryo, embryo);
imp = ij.IJ.openImage(Embryo);
Im = squeeze(ImagePlus2array(imp));

load(sprintf('%s\\Embryo%d\\disp.mat', Folder, embryo));
timepoints = size(disp, 1) + 1;
midline = uint8(timepoints / 2);

orig = struct('x', 0, 'y', 0, 'z', 0);
spacing = struct('x', 1, 'y', 1, 'z', 1);
rifCC = Im(:,:,:);
Seg(:,:,:,midline) = rifCC(:,:,:);

% Creating live-image mask 
tic;
for t = 1:midline-1
    newimCC = mirt3D_transform(double(rifCC), disp{midline-t, 1});
    rifCC = newimCC;
    Seg(:,:,:,midline-t) = newimCC;
end   
toc;

rifCC(:,:,:) = Seg(:,:,:,midline);

tic;
for t = midline:timepoints-1
    newimCC = mirt3D_transform(double(rifCC), disp{t, 1});
    rifCC = newimCC;
    Seg(:,:,:,t+1) = newimCC;
end   
toc;

Seg(Seg(:,:,:,:) == 1) = 255;
imp = copytoImagePlus(Seg, 'YXZT');
ij.IJ.saveAsTiff(imp, sprintf('%s\\Embryo%d\\Segmentation\Seg_completo.tif', Folder, embryo));


%% 2.2.4 Import the segmentation in Fiji. Fill gap manually + close holes --> Fill_(embryo).tif.6 Secting only staged frame from the live-shape mask
clear; close all; clc;

Folder = '\\1. EstimatingIndividualLiveImageMotion';
Excel = xlsread('\\2. IntegratingMultipleLiveImagesIntoAConsensusTemporalReference\EmbryoStage.xlsx'); 

for cla = 1:size(Excel, 1)
    % Load file 
    mask = sprintf('%s\\Embryo%d\\Segmentation\\Seg_completo.tif', Folder, Excel(cla, 2));
    imp = ij.IJ.openImage(mask);
    Im = squeeze(ImagePlus2array(imp));

    time1 = Excel(cla, 3);
    time2 = Excel(cla, 4);
    gr = Excel(cla, 1);
    gr_1 = gr - 1;
    options.overwrite = true;

    Folder1 = sprintf('\\2. IntegratingMultipleLiveImagesIntoAConsensusTemporalReference\\Cluster%d\\Images\\Embryo\\e%d.tif', gr, Excel(cla, 2));
    Folder2 = sprintf('\\2. IntegratingMultipleLiveImagesIntoAConsensusTemporalReference\\Cluster%d\\Images\\Embryo\\e%d.tif', gr_1, Excel(cla, 2));

    saveastiff(Im(:,:,:,time1), Folder1, options);
    saveastiff(Im(:,:,:,time2), Folder2, options);
end



%% 2.2.7 Non-rigid Registration with MIRT. Trasform live-shape mask in ATLAS_Cut mask
clear; close all; clc;
Excel = xlsread('\\2. IntegratingMultipleLiveImagesIntoAConsensusTemporalReference\EmbryoStage.xlsx'); 

for cla = 1:size(Excel, 1)
    Folder = sprintf('\\2. IntegratingMultipleLiveImagesIntoAConsensusTemporalReference\\Cluster%d', Excel(cla, 1));
    embryo = Excel(cla, 2);
    time = Excel(cla, 4);

    SegTest = loadtiff(sprintf('%s\\Images\\Embryo\\%d.%d.tif', Folder, embryo, time));
    SegRef = loadtiff(sprintf('%s\\Images\\Atlas\\Fill_%d.tif', Folder, embryo));

    % Main settings
    main.similarity = 'SSD';
    main.subdivide = 4;
    main.okno = 4;
    main.lambda = 0.1;
    main.single = 1;

    % Optimization settings
    optim.maxsteps = 100;
    optim.fundif = 1e-8;
    optim.gamma = 1;
    optim.anneal = 0.9;

    im = mat2gray(SegTest); % Normalize data 0-1
    refim = mat2gray(SegRef);

    [res, newim] = mirt3D_register(refim, im, main, optim); % Apply MIRT algorithm for non-rigid registration 

    save(sprintf('%s\\Disp\\Non-Rigid\\res%d.mat', Folder, embryo), 'res');
    options.overwrite = true;
    saveastiff(newim, sprintf('%s\\Disp\\DeformedEmbryo\\%d.tif', Folder, embryo), options);

    clearvars -except Excel;

    Folder = sprintf('\\2. IntegratingMultipleLiveImagesIntoAConsensusTemporalReference\\Cluster%d', Excel(cla, 1)-1);
    embryo = Excel(cla, 2);
    time = Excel(cla, 3);

    SegTest = loadtiff(sprintf('%s\\Images\\Embryo\\%d.%d.tif', Folder, embryo, time));
    SegRef = loadtiff(sprintf('%s\\Images\\Atlas\\Fill_%d.tif', Folder, embryo));

    % Main settings
    main.similarity = 'SSD';
    main.subdivide = 4;
    main.okno = 4;
    main.lambda = 0.1;
    main.single = 1;

    % Optimization settings
    optim.maxsteps = 100;
    optim.fundif = 1e-8;
    optim.gamma = 1;
    optim.anneal = 0.9;

    im = mat2gray(SegTest); % Normalize data 0-1
    refim = mat2gray(SegRef);

    [res, newim] = mirt3D_register(refim, im, main, optim); % Apply MIRT algorithm for non-rigid registration 

    save(sprintf('%s\\Disp\\Non-Rigid\\res%d.mat', Folder, embryo), 'res');
    options.overwrite = true;
    saveastiff(newim, sprintf('%s\\Disp\\DeformedEmbryo\\%d.tif', Folder, embryo), options);

    clearvars -except Excel;
end



%% 2.2.8 PointCloud Registration: Force point cloud on ATLAS mask edge --> MAPPED folder (SurfaceMap)
clear; close all; clc;

% Load Excel data
Excel = xlsread('\\2. IntegratingMultipleLiveImagesIntoAConsensusTemporalReference\EmbryoStage.xlsx'); 

for cla = 1:size(Excel, 1)
    % Define folder paths and load required files
    clusterNum = Excel(cla, 1);
    embryo = Excel(cla, 2);
    time = Excel(cla, 4);
    Folder = fullfile('\\2. IntegratingMultipleLiveImagesIntoAConsensusTemporalReference', ['Cluster', num2str(clusterNum)]);
    
    % Load point cloud, transformation, embryo mask, and ATLAS mask
    [node, face] = read_ply(fullfile('\\1. EstimatingIndividualLiveImageMotion', ...
        ['Embryo', num2str(embryo), '\Shapes\CC', num2str(time), '.ply']));
    load(fullfile(Folder, 'Disp\Non-Rigid', ['res', num2str(embryo), '.mat']));
    refim = loadtiff(fullfile(Folder, 'Images\Atlas', ['Fill_', num2str(embryo), '.tif']));
    im = loadtiff(fullfile(Folder, 'Images\Embryo', [num2str(embryo), '.', num2str(time), '.tif']));
    
    % Initialize variables
    point1 = zeros(size(node, 1), 3);
    newcells = zeros(size(refim));  % Create an image for transformed point clouds

    % Enforce all points to be within the image boundaries
    node(:, 2) = max(min(node(:, 2), size(im, 1)), 1);
    node(:, 1) = max(min(node(:, 1), size(im, 2)), 1);
    node(:, 3) = max(min(node(:, 3), size(im, 3)), 1);

    % Decompose the grid transformation into its components
    F = mirt3D_F(res.okno); 
    [Xx, Xy, Xz] = mirt3D_nodes2grid(res.X, F, res.okno);
    X = gpuArray(Xx);
    Y = gpuArray(Xy);
    Z = gpuArray(Xz);

    tic
    % Process each point in the node
    for i = 1:size(node, 1)
        cell = zeros(size(im, 1), size(im, 2));
        cell(int16(node(i, 2)), int16(node(i, 1))) = 1; 

        [r, c] = ind2sub(size(cell), find(cell ~= 0));
        [r1, c1] = ind2sub(size(im(:, :, int16(node(i, 3)))), find(im(:, :, int16(node(i, 3))) ~= 0));

        [Idx, ~] = knnsearch([r1, c1], [r, c]);   
        cells = zeros(size(im));
        cells(int16(r1(Idx)), int16(c1(Idx)), int16(node(i, 3))) = 1;

        rim = interp3(double(cells), X, Y, Z, 'linear', NaN);
        rim = gather(rim);
        rim(isnan(rim)) = 0;

        [r, c, v] = ind2sub(size(rim), find(rim ~= 0));
        [r1, c1, v1] = ind2sub(size(refim), find(refim ~= 0));
        [Idx, D] = knnsearch([r1, c1, v1], [r, c, v]);   
        pos = find(D == min(D), 1);

        if ~isempty(pos)
            newcells(int16(r1(Idx(pos))), int16(c1(Idx(pos))), int16(v1(Idx(pos)))) = 1;
            point1(i, :) = [r1(Idx(pos)), c1(Idx(pos)), v1(Idx(pos))];
        end
    end
    toc

    % Post-processing to smooth and remove mismatches
    [x, y] = find(point1 == 0);
    node2 = node(x, :);
    node3 = node;
    node3(x, :) = 0;

    [Idx, ~] = knnsearch(node3, node2);
    point1(x, :) = point1(Idx, :);
    node = point1;

    % Write output to PLY files
    write_ply(node, face, fullfile(Folder, '3D\Mapped', [num2str(embryo), '.ply']));
    [node, face] = read_ply(fullfile(Folder, '3D\Mapped', [num2str(embryo), '.ply']));
    [newnode2, newface2] = surfreorient(node, face);
    write_ply(newnode2, newface2, fullfile(Folder, '3D\Mapped', [num2str(embryo), '.ply']));

    [node1, face1] = read_ply(fullfile(Folder, '3D\Mapped', [num2str(embryo), '.ply']));
    conn = meshconn(face(:, 1:3), size(node, 1));

    n1{1} = node;
    niter = 4;
    for iter = 2:niter+1
        n1{iter} = smoothsurf(n1{iter-1}, [], conn, 1, 0.9, 'laplacian');
    end
    node = n1{end};
    write_ply(node, face, fullfile(Folder, '3D\Mapped', [num2str(embryo), '.ply']));

    clearvars -except Excel
end

%% 2.2.9 Face-to-Face Matching Live-shape vs. ATLAS_Cut
clear; close all; clc;

% Load Excel data
Excel = xlsread('\\2. IntegratingMultipleLiveImagesIntoAConsensusTemporalReference\EmbryoStage.xlsx'); 

for cla = 1:size(Excel, 1)
    clusterNum = Excel(cla, 1);
    embryo = Excel(cla, 2);
    Folder = fullfile('\\1. EstimatingIndividualLiveImageMotion', ['Cluster', num2str(clusterNum)]);

    % Load Atlas data
    [node, face] = read_ply(fullfile(Folder, '3D\Atlas', ['Cut_', num2str(embryo), '.ply'])); % Cut
    [node1, face1] = read_ply(fullfile(Folder, '3D\Atlas', [num2str(embryo), '.ply'])); % Complete

    % Find the cut part
    centroid = meshcentroid(node, face);
    centroid1 = meshcentroid(node1, face1);
    [IdxCUT, ~] = knnsearch(centroid1, centroid);   

    color = zeros(size(face1, 1), 1);
    color(IdxCUT, :) = 1;

    % Find matching between cut and mapped
    [node, face] = read_ply(fullfile(Folder, '3D\Mapped', [num2str(embryo), '.ply'])); % Mapped 
    [node1, face1] = read_ply(fullfile(Folder, '3D\Atlas', [num2str(embryo), '.ply'])); % Complete 

    % Find the matching part
    centroid = meshcentroid(node, face);
    centroid1 = meshcentroid(node1, face1);
    [IdxMatch, D] = knnsearch(centroid, centroid1(IdxCUT, :));   

    % Save matching indices
    save(fullfile(Folder, 'VertexCorrespondence', ['Embryo', num2str(embryo), '\IdxMatch.mat']), 'IdxMatch');
    save(fullfile(Folder, 'VertexCorrespondence', ['Embryo', num2str(embryo), '\IdxCUT.mat']), 'IdxCUT');

    clearvars -except Excel
end
