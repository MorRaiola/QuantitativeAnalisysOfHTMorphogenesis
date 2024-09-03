%% 3.2 Plotting Deformation onto the ATLAS Shape using Face-to-Face Matching
clear; close all; clc;

% Read data from the Excel file
Excel = xlsread('\\tierra.cnic.es\SC\LAB_MT\RESULTADOS\Morena\EmbryoStage.xlsx'); 

% Loop over each row in the Excel file
for cla = 1:size(Excel, 1)
    
    % Define folder paths
    FolderPre = '\\tierra.cnic.es\SC\LAB_MT\RESULTADOS\Morena\';
    Folder = fullfile(FolderPre, 'Mapping', ['Cluster', num2str(Excel(cla, 1))], 'Deformation');
    FolderTierra = fullfile(FolderPre, 'Defor', ['Gr', num2str(Excel(cla, 1))], '');
    DataCC = fullfile(FolderPre, 'Mapping', 'Atlas', ['remesh', num2str(Excel(cla, 1)), '.ply']);
    
    % Load the PLY file containing node and face data
    [node1, face1] = read_ply(DataCC);
    node1 = node1(:, [3, 1, 2]); % Reorder the vertices
    face1 = face1(:, [3, 1, 2]); % Reorder the face indices
    
    % Create a struct for the mesh data
    myo.vertices = node1;
    myo.faces = face1;
    
    % Get the current embryo and cluster numbers
    embryo = Excel(cla, 2);
    cluster = num2str(Excel(cla, 1));
    
    %% MAP GROWTH RATE
    % Load deformation and face-matched indices
    load(fullfile(Folder, ['JCCsmooth', num2str(embryo), '.mat']));
    load(fullfile(FolderPre, 'Mapping', ['Cluster', cluster], 'VertexCorrespondence', ['Embryo', num2str(embryo)], 'IdxCUT.mat'));
    load(fullfile(FolderPre, 'Mapping', ['Cluster', cluster], 'VertexCorrespondence', ['Embryo', num2str(embryo)], 'IdxMatch.mat'));

    JCC = newJCC; % Update the variable for clarity
    color = zeros(size(face1, 1), size(JCC, 2));
    
    % Apply deformation values to the corresponding vertices
    for i = 1:size(color, 2)
        color(IdxCUT, i) = JCC(IdxMatch, i);
    end
    
    % Plot the mesh with growth rate coloring
    plotmesh(node1, face1, 'FaceVertexCData', color, 'EdgeColor', 'none', 'FaceLighting', 'gouraud', ...
        'AmbientStrength', 0.3, 'DiffuseStrength', 0.8, 'SpecularStrength', 0.9, 'SpecularExponent', 25, 'BackFaceLighting', 'unlit');
    view([-180, 90])
    caxis([Excel(cla, 15), Excel(cla, 16)])
    colormap('jet')
    axis equal off
    lightangle(0, 180)
    camzoom(1.5)
    close 
    
    % Delete existing files and save the new OBJ file
    delete(fullfile(Folder, 'Atlas', 'JCC', [num2str(Excel(cla, 2)), '.obj']));
    delete(fullfile(Folder, 'Atlas', 'JCC', [num2str(Excel(cla, 2)), '.mtl']));
    obj_write_color(myo, fullfile(Folder, 'Atlas', 'JCC', num2str(Excel(cla, 2))), color, 'cmin', Excel(cla, 7));
    
    %% MAP ANISOTROPY
    % Load deformation and face-matched indices for anisotropy
    load(fullfile(Folder, ['Anisotropysmooth', num2str(embryo), '.mat']));
    load(fullfile(FolderPre, 'Mapping', ['Cluster', cluster], 'VertexCorrespondence', ['Embryo', num2str(embryo)], 'IdxCUT.mat'));
    load(fullfile(FolderPre, 'Mapping', ['Cluster', cluster], 'VertexCorrespondence', ['Embryo', num2str(embryo)], 'IdxMatch.mat'));

    Anisotropy = newAnisotropy;
    color = zeros(size(face1, 1), size(Anisotropy, 2));
    
    % Apply anisotropy values to the corresponding vertices
    for i = 1:size(color, 2)
        color(IdxCUT, i) = Anisotropy(IdxMatch, i);
    end
    
    % Plot the mesh with anisotropy coloring
    plotmesh(node1, face1, 'FaceVertexCData', color, 'EdgeColor', 'none', 'FaceLighting', 'gouraud', ...
        'AmbientStrength', 0.3, 'DiffuseStrength', 0.8, 'SpecularStrength', 0.9, 'SpecularExponent', 25, 'BackFaceLighting', 'unlit');
    view([-180, 90])
    colormap('jet')
    caxis([Excel(cla, 17), Excel(cla, 18)])
    axis equal off
    lightangle(0, 180)
    camzoom(1.5)
    
    % Delete existing files and save the new OBJ file
    delete(fullfile(Folder, 'Atlas', 'Anisotropy', [num2str(Excel(cla, 2)), '.obj']));
    delete(fullfile(Folder, 'Atlas', 'Anisotropy', [num2str(Excel(cla, 2)), '.mtl']));
    obj_write_color(myo, fullfile(Folder, 'Atlas', 'Anisotropy', num2str(Excel(cla, 2))), color, 'cmin', Excel(cla, 9));

    %% MAP STRAIN
    % Load deformation and face-matched indices for strain
    load(fullfile(Folder, ['StrainMaxsmooth', num2str(embryo), '.mat']));
    load(fullfile(FolderPre, 'Mapping', ['Cluster', cluster], 'VertexCorrespondence', ['Embryo', num2str(embryo)], 'IdxCUT.mat'));
    load(fullfile(FolderPre, 'Mapping', ['Cluster', cluster], 'VertexCorrespondence', ['Embryo', num2str(embryo)], 'IdxMatch.mat'));

    Strain = newStrain;
    color = zeros(size(face1, 1), size(Strain, 2));
    
    % Apply strain values to the corresponding vertices
    for i = 1:size(color, 2)
        color(IdxCUT, i) = Strain(IdxMatch, i);
    end
    
    % Save the strain data
    save(fullfile(Folder, 'Atlas', 'Strain', ['Stretch', num2str(embryo), '.mat']), 'color');
    
    % Plot the mesh with strain coloring
    plotmesh(node1, face1, 'FaceVertexCData', color, 'EdgeColor', 'none', 'FaceLighting', 'gouraud', ...
        'AmbientStrength', 0.3, 'DiffuseStrength', 0.8, 'SpecularStrength', 0.9, 'SpecularExponent', 25, 'BackFaceLighting', 'unlit');
    view([-180, 90])
    colormap('jet')
    caxis([Excel(cla, 19), Excel(cla, 20)])
    axis equal off
    lightangle(0, 180)
    camzoom(1.5)
    
    % Delete existing files and save the new OBJ file
    delete(fullfile(Folder, 'Atlas', 'Strain', [num2str(Excel(cla, 2)), '.obj']));
    delete(fullfile(Folder, 'Atlas', 'Strain', [num2str(Excel(cla, 2)), '.mtl']));
    obj_write_color(myo, fullfile(Folder, 'Atlas', 'Strain', num2str(Excel(cla, 2))), color, 'cmin', Excel(cla, 11));
    
    % Clear variables except for Excel and loop variables
    clearvars -except Excel valueStr valueAni valueJCC
end
