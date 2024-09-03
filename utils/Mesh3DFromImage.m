function [newnode2, newface2]=Mesh3DFromImage(inputTiffPath, outputPlyPath, plotFigures, saveimage)
    % Mesh3DFromImage - Generates a 3D mesh from a segmented TIFF image and saves it as a PLY file.
    %
    % Syntax:
    %    Mesh3DFromImage(inputTiffPath, outputPlyPath, plotFigures)
    %
    % Inputs:
    %    inputTiffPath - Full path to the input segmented TIFF file.
    %    outputPlyPath - Full path where the output PLY file will be saved.
    %    plotFigures   - Optional boolean flag to enable/disable plotting (default is true).
    %
    % Example:
    %    Mesh3DFromImage('G:\Embryos\Embryo1\Segmentation\Embryo1.tif', ...
    %                    'C:\Users\mraiola\Desktop\Endo\SegResize.ply', true);

    if nargin < 3
        plotFigures = true; % Default to plotting if not specified
    end

    if nargin < 4
        plotFigures = true; % Default to plotting if not specified
    end


    % Load the TIFF segmentation file
    Im = loadtiff(inputTiffPath);

    % Isolate myocardium by retaining only label 1 and zeroing others
    Im(Im == 1) = 1;  % Retain myocardium
    Im(Im ~= 1) = 0;  % Set all other segments to zero

    % Fill holes in the segmented 3D volume and generate the mesh
    cleancc = fillholes3d(logical(Im > 0), 1); 
    [node, ~, face] = v2m(cleancc, 0.1, 5, 5);
   
    % Plot the initial mesh if plotting is enabled
    if plotFigures
        figure;
        plotmesh(node, face, 'facealpha', 0.7);
        title('Initial Mesh');
    end

    % Smooth the surface using a Laplacian filter
    conn = meshconn(face(:, 1:3), size(node, 1));
    n1{1} = node;
    niter = 6;  % Number of smoothing iterations

    for iter = 2:niter + 1
        n1{iter} = smoothsurf(n1{iter - 1}, [], conn, 1, 0.9, 'laplacian');
    end

    % Reorient and finalize the smoothed mesh
    node = n1{end}(:, [2, 1, 3]);  % Reorder coordinates
    [newnode2, newface2] = surfreorient(node, face); 

    % Plot the smoothed mesh if plotting is enabled
    if plotFigures
        figure;
        plotmesh(newnode2, newface2, 'facealpha', 0.7);
        title('Smoothed Mesh');
    end 

    % Save the final mesh to the specified PLY file
    if saveimage
        write_ply(newnode2, newface2, outputPlyPath);
    end
end
