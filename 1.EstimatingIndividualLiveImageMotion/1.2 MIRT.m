%% 1.2 Estimating Individual Live Image Motion using Mirt3D Algorithm at t(N/2)
clear; close all; clc;

% Get the folder path of the current script
folder = fileparts(mfilename('fullpath'));

% Add MIJ and ImageJ Java paths
addpath('C:\Users\mraiola\Downloads\fiji-win64\Fiji.app\scripts');
javaaddpath('C:\Program Files\MATLAB\R2020a\java\jar\mij.jar');
javaaddpath('C:\Program Files\MATLAB\R2020a\java\jar\ij-1.52i.jar');

% Load the time-lapse image data
Embryo = '\1. Estimating Individual Live Image Motion\Embryo1\Data\Embryo1.tif';
Output = '\1. Estimating Individual Live Image Motion\Embryo1';

imp = ij.IJ.openImage(Embryo); 
Im = ImagePlus2array(imp);
Im = squeeze(Im); % Reshape image to YXZT format

% Main settings for the Mirt3D algorithm
main.similarity = 'SSD';   % Similarity measure (Sum of Squared Differences)
main.subdivide  = 4;       % Number of hierarchical levels
main.okno       = 4;       % Mesh grid resolution
main.lambda     = 0.01;    % Regularization weight for transformation
main.single     = 0;       % Display mesh transformation at each iteration (set to 1 for debugging)

% Optimization settings for the Mirt3D algorithm
optim.maxsteps = 100;      % Maximum number of iterations at each level
optim.fundif   = 1e-8;     % Tolerance for stopping criterion
optim.gamma    = 1;        % Initial optimization step size
optim.anneal   = 0.9;      % Annealing rate for step size reduction

% Initialize displacement fields
orig.x = 0; spacing.x = 1;
orig.y = 0; spacing.y = 1;
orig.z = 0; spacing.z = 1;

% Set up parallel processing
parpool(6);

% Determine the number of timepoints and the midpoint (t = N/2)
timepoints = size(Im, 4);
midline = uint8(timepoints / 2);


% Forward pass: Register images from the midpoint to the first frame
tic
parfor t = 1:midline-1
    im = mat2gray(Im(:,:,:,midline - t + 1));
    refim = mat2gray(Im(:,:,:,midline - t));
    [res, ~] = mirt3D_register(refim, im, main, optim);
    disp{t,1} = res;
end

% Reverse the order of displacement fields for the first half
disp2 = cell(timepoints - 1, 1);
for k = 1:size(disp, 1)
    disp2{midline - k, 1} = disp{k, 1};
end

% Backward pass: Register images from the midpoint to the last frame
parfor t = midline:timepoints - 1
    im = mat2gray(Im(:,:,:,t));
    refim = mat2gray(Im(:,:,:,t + 1));   
    [res, ~] = mirt3D_register(refim, im, main, optim);
    disp2{t,1} = res;
end
toc

% Save the displacement field results
save([Output 'transformation.mat'], 'disp2');
