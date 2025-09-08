%% Modelling the Staging System
clear; close all; clc

% Define the filename and read data from the Excel file
filename = 'V:\PaperFigureData\Figure3_Integrating_Multiple_Live_Images_Into_A_Consensus_Spatio_Temporal_Reference\Figure3.xlsx';
GNB = readmatrix(filename, 'Sheet', 'B');  % Read data from the sheet named 'B'
liveimages = readmatrix(filename, 'Sheet', 'C');  % Read data from the sheet named 'C'

% Extract Gaussian parameters from GNB data
gr = GNB(2:10, 7);    % Group labels or identifiers
mu = GNB(2:10, 8);    % Means of Gaussian distributions
sigma = GNB(2:10, 9); % Standard deviations of Gaussian distributions

% Extract features from liveimages and round them
features = liveimages(:, 1:5:end);
liveimages(:, 1:5:end) = round(features, 5);

% Define column indices for group and probability values in liveimages
findGr = 4:5:size(liveimages, 2) + 1;
findPP = 5:5:size(liveimages, 2) + 2;

% Process each feature set
for e = 1:size(features, 2)
    singleemb = features(:, e);        % Extract single feature set
    singleemb(isnan(singleemb)) = [];  % Remove NaN values

    % Compute likelihoods using the Gaussian PDFs
    y = normpdf(singleemb, mu', sigma');

    % Normalize the probabilities
    S = sum(y, 2);  % Compute the sum of all Gaussian PDFs for normalization
    table = y ./ S; % Normalize PDFs to obtain probabilities

    % Find the maximum probability and corresponding Gaussian index for each data point
    [M, I] = max(table, [], 2);
    I = gr(I);  % Map indices to group labels

    % Update the liveimages matrix with the rounded group labels and probabilities
    liveimages(1:size(I, 1), findGr(1, e)) = round(I, 4);
    liveimages(1:size(M, 1), findPP(1, e)) = round(M, 4);

    % Clear temporary variables
    clear y table M I 
end

% Define the output filename and write the results to an Excel file
outputFilename = 'F:\clustering1.xlsx';
writematrix(liveimages, outputFilename);
