%% 1.1 Processing Embryo Image Data
clear; close all; clc;

% Get the folder path of the current script
folder = fileparts(mfilename('fullpath')); 

% Add MIJ and ImageJ Java paths
addpath('C:\Users\mraiola\Downloads\fiji-win64\Fiji.app\scripts');
javaaddpath('C:\Program Files\MATLAB\R2020a\java\jar\mij.jar'); % MIJ Java folder
javaaddpath('C:\Program Files\MATLAB\R2020a\java\jar\ij-1.52i.jar'); % ImageJ JAR file

% Define input and output paths
Embryo = '... \1. Estimating Individual Live Image Motion\Embryo1\Data\Embryo1.tif'; % Path to the live image file
Output = '... \1. Estimating Individual Live Image Motion\Embryo1\Data\'; % Output directory

% Load the image using ImageJ
imp = ij.IJ.openImage(Embryo);
Im = ImagePlus2array(imp);
Im = squeeze(Im); % Reshape image to YXZT format

% Define pre-processing parameters
ImRes = [0.593, 0.593, 5]; % Image resolution in micrometers (µm)
downS = 0.25; % Downscaling factor
crop = [54, 114, 1, 248, 170, 137]; % Cuboid region to crop: [x, y, z, width, height, depth]

% Initialize an array for the processed image
EmbryoResize = [];

% Process each time point (4th dimension in the image stack)
for t = 1:size(Im, 4)
    % Crop the image to the specified region
    EmbryoCropped(:,:,:,t) = imcrop3(Im(:,:,:,t), crop);
    
    % Move data to the GPU for filtering (optional)
    im = gpuArray(EmbryoCropped(:,:,:,t)); 
    
    % Apply Gaussian filtering to smooth the image
    EmbryoGauss = imgaussfilt3(im, 1);
    
    % Move the data back to the CPU
    outCPU = gather(EmbryoGauss);
    
    % Apply median filtering to remove noise
    EmbryoMedian = medfilt3(outCPU);
    
    % Resize the image based on the downscaling factor and resolution ratio
    EmbryoResize(:,:,:,t) = imresize3(EmbryoMedian, ...
        [size(EmbryoMedian,1)*downS, size(EmbryoMedian,2)*downS, size(EmbryoMedian,3)*ImRes(3)/ImRes(1)*downS]);
end 

% Convert the processed image back to an ImagePlus object for saving
imp = copytoImagePlus(EmbryoResize, 'YXZT');

% Save the processed image as a TIFF file
ij.IJ.saveAsTiff(imp, fullfile(Output, 'Embryo1.tif'));

