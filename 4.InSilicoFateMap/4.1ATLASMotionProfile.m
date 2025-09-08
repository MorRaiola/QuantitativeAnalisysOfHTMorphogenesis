%% Concatenation to build the Fate Map
%% 4.1 Continous Model 
clear; close all; clc

% Define cut and embryo arrays
cutArray = { [3, 4, 5, 6, 7], 2, [8, 9] }; % Array of cuts
embryoArray = { [31, 16, 35, 27, 24, 12], 31, 12 }; % Corresponding embryos
clusterFolder = '\2. IntegratingMultipleLiveImagesIntoAConsensusTemporalReference\Cluster';
deformationFolder = '\3. QuantifyingTissueDeformation\DinamicAtlas\';

% Loop through each set of cuts and embryos
for i = 1:numel(cutArray)
    cuts = cutArray{i};
    embryos = embryoArray{i};
    
    for cutIndex = 1:numel(cuts)
        % Adjust embryo index based on the current cut
        embryoIndex = min(cutIndex + (i - 1) * 2, numel(embryos));
        
        % Read PLY file and .mat file
        [node, face] = read_ply([clusterFolder num2str(cuts(cutIndex)) ...
                                   '\3D\Mapped\' num2str(embryos(embryoIndex)) '.ply']);
        load([clusterFolder num2str(cuts(cutIndex)) ...
              '\Disp\Rigid\' num2str(embryos(embryoIndex)) '.mat']);
        
        % Apply transformation to nodes
        tempSR = bsxfun(@times,Transform.R,Transform.s);
        trPts = (tempSR*node');
        Xo = bsxfun(@plus,trPts',Transform.t);
                
        write_ply(transformedPoints', face, [deformationFolder ...
            num2str(cuts(cutIndex)) '\' num2str(embryos(embryoIndex)) 'Rot.ply']);
    end
end

%% 4.2 Cut (embryo)Rot.ply according e31.ply, since e31 represents our common tissue part. Result -> (embryo)Rot_Cut.ply


%% 4.3 Since index are not respected during Meshlab cutting we replaced them 
clear;close all;clc

cutArray = { [3, 4, 5, 6, 7], 2, [8, 9] }; % Array of cuts
embryoArray = { [31, 16, 35, 27, 24, 12], 31, 12}; % Corresponding embryos
clusterFolder = '\2. IntegratingMultipleLiveImagesIntoAConsensusTemporalReference\Cluster';
deformationFolder = '\3. QuantifyingTissueDeformation\DinamicAtlas\';

for i = 1:numel(cutArray)
    cuts = cutArray{i};
    embryos = embryoArray{i};
    
    for cutIndex = 1:numel(cuts)
                embryoIndex = min(cutIndex + (i - 1) * 2, numel(embryos));

[node,face] = read_ply([deformationFolder num2str(cuts(cutIndex)) '\' num2str(embryos(embryoIndex)) 'Rot_Cut.ply']); %% meshlab results
[node1,face1] = read_ply([clusterFolder num2str(cuts(cutIndex)) '\3D\Mapped\' num2str(embryos(embryoIndex)) '.ply']);
load([clusterFolder num2str(cuts(cutIndex)) '\Disp\Rigid\' num2str(embryos(embryoIndex)) '.mat']);
tempSR = bsxfun(@times,Transform.R,Transform.s);
trPts = (tempSR*node1');
Xo = bsxfun(@plus,trPts',Transform.t);
centroid = meshcentroid(node,face);
centroid1 = meshcentroid(Xo,face1);
[IdxMatch,D] = knnsearch(centroid1,centroid);  
face1 = face1(IdxMatch,:);
write_ply(Xo,face1,[deformationFolder num2str(cuts(cutIndex)) '\' num2str(embryos(embryoIndex)) 'Rot_Cut.ply']); %% with the right index
    end
end


%% 4.4 Describing tissue with the same num of points (we registered with 35's num of points)
clear; close all; clc;
deformationFolder = '\3. QuantifyingTissueDeformation\DinamicAtlas\';
clusterFolder = '\2. IntegratingMultipleLiveImagesIntoAConsensusTemporalReference\Cluster';

% Define the cuts and corresponding embryos and times for processing
data = {
    [4, 5], [16, 27], [3, 4, 0, 6]; 
    3, 31, [2, 3];                 
    6, 24, [6, 7];                
    7, 12, [7, 8, 9]              
};

% Loop through each dataset
for d = 1:size(data, 1)
    cut = data{d, 1}; 
    embryo = data{d, 2}; 
    time = data{d, 3}; 
    pos = find(time == 0); 
    s = 1; m = time(1);
    
    for c = 1:numel(cut)
        [node, face] = read_ply([deformationFolder  num2str(cut(c)) '\35Rot_Cut.ply']);
        [node1, face1] = read_ply([deformationFolder  num2str(cut(c)) '\' num2str(embryo(1, c)) 'Rot_Cut.ply']);
        
        % Set registration options
        opt.method='rigid'; % use rigid registration
        opt.viz=1;          % show every iteration
        opt.outliers=0.6;   % use 0.6 noise weight to add robustness 
        
        opt.normalize=1;    % normalize to unit variance and zero mean before registering (default)
        opt.scale=1;        % estimate global scaling too (default)
        opt.rot=1;          % estimate strictly rotational matrix (default)
        opt.corresp=1;      % do not compute the correspondence vector at the end of registration (default)
        
        opt.max_it=100;     % max number of iterations
        opt.tol=1e-8;       % tolerance

        % Register Y to X
        [Transform, Correspondence] = cpd_register(node1, node, opt);
        write_ply(node1(Correspondence,:), face, [deformationFolder num2str(cut(c)) '\' num2str(embryo(1, s)) 'Rot_Registered.ply']);
        
        if c==size(pos,2)+1
            subtime = m:time(end);
        else
            subtime = m:time(pos(c)-1);
        end

        for k = 1:numel(subtime)
            [node1, face1] = read_ply([clusterFolder num2str(subtime(k)) '\3D\Mapped\' num2str(embryo(1, s)) '.ply']);
            write_ply(node1(Correspondence,:), face, [deformationFolder num2str(subtime(k)) '\' num2str(embryo(1, s)) 'Rot_Registered.ply']);
        end
        
        s = s + 1; m = time(pos + 1);
    end
end

%% 4.5 To map the in silico fate Map into the ATLAS_ Defining the IdxCut and IdxMatch betwent fate map and ATLAS
clear; close all; clc

% Loop through classes 2 to 9
for cla = 2:9
    Folder = ['\\3. QuantifyingTissueDeformation\DinamicAtlas\' num2str(cla)];
    inFiles = dir([Folder filesep '*.ply']); 
    inNames = {};
    
    % Collect all .ply file names in the folder
    for i = 1:length(inFiles)
        inNames{end+1} = [Folder filesep inFiles(i).name];
    end
    
    % Extract numerical identifiers from file names
    if size(inNames,2)>8
       [~,name,~] = fileparts(inNames{1});
       num(1,1) = str2double(name);
       [filepath,name,ext] = fileparts(inNames{6});
       num(2,1) = str2double(name);
    else
       [filepath,name,ext] = fileparts(inNames{1});
       num(1,1) = str2double(name);   
    end

    % Load atlas complete and atlas cut
    for i = 1:size(num, 1)
        % Read cut and complete meshes
        [node, face] = read_ply([Folder '\' num2str(num(i)) '_Cut.ply']);   % Cut mesh
        [node1, face1] = read_ply([Folder '\' num2str(num(i)) '.ply']);    % Complete mesh 

        % Find the centroid of the meshes
        centroid = meshcentroid(node, face);
        centroid1 = meshcentroid(node1, face1);
        [IdxCUT, ~] = knnsearch(centroid1, centroid);    % Find nearest neighbors

        % Create color array to visualize the cut part
        color = zeros(size(face1, 1), 1);
        color(IdxCUT) = 1;  % Mark cut part

        % Uncomment to visualize the cut part
        % figure
        % plotmesh(node1, face1, 'FaceVertexCData', color, 'FaceColor', 'flat', ...
        %          'EdgeColor', 'none', 'FaceLighting', 'gouraud', 'AmbientStrength', 0.5);

        % Find the matching between cut and mapped meshes
        [node, face] = read_ply([Folder '\' num2str(num(i)) 'Rot_Registered.ply']); % Mapped mesh
        [node1, face1] = read_ply([Folder '\' num2str(num(i)) '.ply']);             % Complete mesh 

        % Register using ICP
        [tform, movingReg] = pcregistericp(pointCloud(node), pointCloud(node1));

        % Find the matching part
        centroid = meshcentroid(movingReg.Location, face);
        centroid1 = meshcentroid(node1, face1);
        [IdxMatch, D] = knnsearch(centroid, centroid1(IdxCUT, :));   

        % Save indices of matches
        save([Folder '\' num2str(num(i)) 'IdxMatch.mat'], 'IdxMatch');
        save([Folder '\' num2str(num(i)) 'IdxCUT.mat'], 'IdxCUT');
    end
end



