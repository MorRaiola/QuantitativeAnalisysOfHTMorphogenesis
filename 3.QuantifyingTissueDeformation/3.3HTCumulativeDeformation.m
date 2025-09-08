%% Concatenating path. 
% we set as reference mesh the e27 SurfaceMaps

%1 step, crea immagini della SurfaceMap.
clear; close all; clc

% Define parameters for each group
group_data = {
    'Cluster4', [2 7 16 22 27 35 124], 'Fill_27.tif';
    'Cluster3', [5 22 31], 'Fill_22.tif';
    'Cluster6', [24 27], 'Fill_27.tif';
    'Cluster8', [12 24], 'Fill_24.tif'
};

baseFolder = '\Code_data\2. IntegratingMultipleLiveImagesIntoAConsensusTemporalReference\';

for group_idx = 1:size(group_data, 1)
    clusterName = group_data{group_idx, 1};
    embryos = group_data{group_idx, 2};
    atlasFile = group_data{group_idx, 3};

    FolderPre = fullfile(baseFolder, clusterName);

    for cla = 1:numel(embryos)
        OutFolder = fullfile(FolderPre, 'Deformation', 'DynamicAtlas', [num2str(embryos(cla)) '.tif']);
        [node, face] = read_ply(fullfile(FolderPre, '3D', 'Mapped', [num2str(embryos(cla)) '.ply']));
        SegTest = loadtiff(fullfile(FolderPre, 'Images', 'Atlas', atlasFile));

        dx = 1:size(SegTest, 1);
        dy = 1:size(SegTest, 2);
        dz = 1:size(SegTest, 3);

        imag = surf2volz(node, face, dx, dy, dz);
        options.overwrite = true;
        saveastiff(imag, OutFolder, options);
    end
    
    % Clear variables for the next group iteration
    clearvars -except group_data baseFolder
end


%% Rigid registration of the image to oriented and resize the image 
clear; close all; clc

% Define parameters for each group (same structure as in previous code)
group_data = {
    'Cluster4', [2 7 16 22 35 124], '27.tif';
    'Cluster3', [5 31], '22.tif';
    'Cluster6', [24], '27.tif';
    'Cluster8', [12], '24.tif'
};

baseFolder = '\Code_data\2. IntegratingMultipleLiveImagesIntoAConsensusTemporalReference';

for group_idx = 3:size(group_data, 1)
    clusterName = group_data{group_idx, 1};
    embryos = group_data{group_idx, 2};
    refFile = group_data{group_idx, 3};

    FolderPre = fullfile(baseFolder, clusterName);

    for e = 1:numel(embryos)
        % Load SegTest and SegRef images
        SegTest = loadtiff(fullfile(FolderPre, 'Deformation', 'DynamicAtlas', [num2str(embryos(e)) '.tif']));
        SegRef = loadtiff(fullfile(FolderPre, 'Deformation', 'DynamicAtlas', refFile));

        % Convert to uint8 if needed
        SegTest = uint8(SegTest);
        SegRef = uint8(SegRef);

        % Display comparison between the center slices
        centerFixed = size(SegRef, 3) / 2;
        centerMoving = size(SegTest, 3) / 2;
        imshowpair(SegTest(:, :, uint8(centerMoving)), SegRef(:, :, uint8(centerFixed)));

        % 3D Image registration
        Rfixed = imref3d(size(SegRef), 1, 1, 1);
        Rmoving = imref3d(size(SegTest), 1, 1, 1);
        [optimizer, metric] = imregconfig('multimodal');
        optimizer.InitialRadius = 0.00063;
        tform = imregtform(SegTest, Rmoving, SegRef, Rfixed, "similarity", optimizer, metric);
        movingRegistered = imwarp(SegTest, Rmoving, tform, "bicubic", 'OutputView', Rfixed);

        % Save registered image and transformation matrix
        options.overwrite = true;
        saveastiff(movingRegistered, fullfile(FolderPre, 'Deformation', 'DynamicAtlas', ['resize' num2str(embryos(e)) '.tif']), options);
        save(fullfile(FolderPre, 'Deformation', 'DynamicAtlas', ['RigidTransf' num2str(embryos(e)) '.mat']), 'tform');

        % Load and apply transformation to 3D surface mesh
        [node, face] = read_ply(fullfile(FolderPre, '3D', 'Mapped', [num2str(embryos(e)) '.ply']));
        ptCloudOut = pctransform(pointCloud(node), tform);
        write_ply(ptCloudOut.Location, face, fullfile(FolderPre, 'Deformation', 'DynamicAtlas', ['resize' num2str(embryos(e)) '.ply']));
    end
    
    % Clear variables for the next group iteration
    clearvars -except group_data baseFolder outputFolder
end


%% Non-rigid registration.
clear; close all; clc

% Define group data (cluster and embryos)
group_data = {
    'Cluster4', [2 7 16 22 35 124], '27.tif';
    'Cluster3', [5 31], '22.tif';
    'Cluster6', [24], '27.tif';
    'Cluster8', [12], '24.tif'
};

baseFolder = '\Code_data\2. IntegratingMultipleLiveImagesIntoAConsensusTemporalReference;

% Main settings
main.similarity = 'SSD';   % similarity measure, e.g. SSD, CC, SAD, RC, CD2, MS, MI
main.subdivide  = 4;       % use 3 hierarchical levels
main.okno       = 4;       % mesh window size
main.lambda     = 0.1;     % transformation regularization weight, 0 for none
main.single     = 1;       % show mesh transformation at every iteration

% Optimization settings
optim.maxsteps = 100;      % maximum number of iterations at each hierarchical level
optim.fundif   = 1e-8;     % tolerance (stopping criterion)
optim.gamma    = 1;        % initial optimization step size 
optim.anneal   = 0.9;      % annealing rate on the optimization step

for group_idx = 1:size(group_data, 1)
    clusterName = group_data{group_idx, 1};
    embryos = group_data{group_idx, 2};
    refFile = group_data{group_idx, 3};

    FolderPre = fullfile(baseFolder, clusterName);
    
    % Initialize results and images
    disp = cell(numel(embryos), 1);  
    images = cell(numel(embryos), 1);  

    % Parallel loop for processing each embryo
    parfor e = 1:numel(embryos)
        % Load the test and reference images
        SegTest = loadtiff(fullfile(FolderPre, 'Deformation', 'DynamicAtlas', ['resize' num2str(embryos(e)) '.tif']));
        SegRef = loadtiff(fullfile(FolderPre, 'Deformation', 'DynamicAtlas', refFile));
        
        % Normalize images to [0,1]
        im = mat2gray(SegTest);
        refim = mat2gray(SegRef);
        
        % Perform 3D registration
        [res, newim] = mirt3D_register(refim, im, main, optim);
        
        % Store results
        disp{e,1} = res;
        images{e,1} = newim;
    end

    % Save the results and registered images
    for e = 1:numel(embryos)
        res = disp{e,1};
        newim = images{e,1};
        
        % Save transformation result
        save(fullfile(FolderPre, 'Deformation', 'DynamicAtlas', [num2str(embryos(e)) '.mat']), 'res');
        
        % Save registered image
        options.overwrite = true;
        saveastiff(newim, fullfile(FolderPre, 'Deformation', 'DynamicAtlas', ['registered' num2str(embryos(e)) '.tif']), options);
    end
    
    % Clear variables before next group
    clearvars -except group_data baseFolder main optim
end


%% Register Groups: Otput = Registered.ply
Folder = '\Code_data\2. IntegratingMultipleLiveImagesIntoAConsensusTemporalReference\';
SubFolder = '\Deformation\DynamicAtlas\';

% Register different groups
registerGroup([2 7 16 22 35 124], '27', 'Cluster4', Folder, SubFolder);
registerGroup([5 31], '22', 'Cluster3', Folder, SubFolder);
registerGroup([24], '27', 'Cluster6', Folder, SubFolder);
registerGroup([12], '24', 'Cluster8', Folder, SubFolder);

%% Face-to-Face matching
clear; close all; clc
embryoSet = {[2 7 16 22 35 124], [24]};
clusterSet = [4, 6];
plyFileSet = {'27.ply', '27.ply'};
FolderPre = '\Code_data\2. IntegratingMultipleLiveImagesIntoAConsensusTemporalReference\Cluster';
SubFolder = '\Deformation\DynamicAtlas\';

for i = 1:length(embryoSet)
    embryos = embryoSet{i};
    cluster = clusterSet(i);
    plyFile = plyFileSet{i};
    
    for e = 1:numel(embryos)
        [nodeRef, faceRef] = read_ply([FolderPre num2str(cluster) SubFolder plyFile]);
        IdxMatch = knnsearch(node, nodeRef);
        save([FolderPre num2str(cluster) SubFolder 'IdxMatch' num2str(embryos(e)) '.mat'], 'IdxMatch');

        centroid = meshcentroid(node, face);
        centroidRef = meshcentroid(nodeRef, faceRef);
        IdxMatchCentroid = knnsearch(centroid, centroidRef);
        save([FolderPre num2str(cluster) SubFolder 'IdxMatchCentroid' num2str(embryos(e)) '.mat'], 'IdxMatchCentroid');
    end
end



%% Propagation Gr4 (e27)
clear;close all; clc
embryo = [2,35,16,124,22];
cluster = 4;
anisotropy = 0;
FolderPre = '\Code_data\2. IntegratingMultipleLiveImagesIntoAConsensusTemporalReference\';
SubFolder = '\Deformation\DynamicAtlas\';

for e = 1:size(embryo,2)
    if anisotropy==1
        name = 'colormeanAnisotropy';
        singlecolor = 'colorAni';
        defor ='Anisotropysmooth';
    else
        name = 'SurfArea';
        singlecolor = 'colorGroSurf';
        defor = 'SurfArea';
    end

[node,face] = read_ply([FolderPre num2str(cluster) SubFolder 'registered' num2str(embryo(e)) '.ply']);
conn = meshconn(face(:,1:3),size(node,1));
n1{1} = node;
niter = 2;

    for iter = 2:niter+1
      n1{iter} = smoothsurf(n1{iter-1},[],conn,1,0.9,'laplacian');
    end 

node = n1{end};

[nodeRef,faceRef] = read_ply([FolderPre num2str(cluster) '\3D\Mapped\27.ply']);
load([FolderPre num2str(cluster) SubFolder '\IdxMatch' num2str(embryo(e)) '.mat']);
load([FolderPre num2str(cluster) '\Deformation\' defor num2str(embryo(e)) '.mat']);
load([FolderPre num2str(cluster) SubFolder '\IdxMatchCentroid' num2str(embryo(e)) '.mat']);


    if anisotropy==1
        color = newAnisotropy(IdxMatchCentroid); 
    else
        color = areas(IdxMatchCentroid); 
    end

myo.vertices = nodeRef;
myo.faces = faceRef;

delete ([FolderPre num2str(cluster) SubFolder singlecolor num2str(embryo(e)) '.mtl'])
delete ([FolderPre num2str(cluster) SubFolder singlecolor num2str(embryo(e)) '.obj'])
obj_write_color(myo, [FolderPre num2str(cluster) SubFolder singlecolor num2str(embryo(e))],color)

colormerge(:,e) = color;

    if e == size(embryo,2)
        colormean = mean(colormerge,2);
        save([FolderPre num2str(cluster) SubFolder name '.mat'],'colormean');
        delete ([FolderPre num2str(cluster) SubFolder name '.mtl'])
        delete ([FolderPre num2str(cluster) SubFolder name '.obj'])
        obj_write_color(myo, [FolderPre num2str(cluster) SubFolder name],colormean)
        write_ply(nodeRef,faceRef,[FolderPre num2str(cluster) SubFolder 'mean.ply']);
    end
end

%% Propagation Gr3 (e22)
clear;close all; clc
embryo= [22,16,124];
cluster = 3;
anisotropy = 0;
FolderPre = '\Code_data\2. IntegratingMultipleLiveImagesIntoAConsensusTemporalReference\Cluster';
SubFolder = '\Deformation\DynamicAtlas\';

for e = 1:size(embryo,2)
    if anisotropy==1
        name = 'colormeanAnisotropy';
        singlecolor = 'colorAni';
        defor = 'Anisotropysmooth';
    else
        name = 'colormean';
        singlecolor = 'colorGro';
        defor ='JCCsmooth';
    end

[node,face] = read_ply([FolderPre num2str(cluster) '\3D\Mapped\' num2str(embryo(e)) '.ply']);
[nodeRef,faceRef] = read_ply([FolderPre num2str(cluster+1) '\3D\Mapped\27.ply']);
load([FolderPre num2str(cluster) '\Deformation\' defor num2str(embryo(e)) '.mat']);
load([FolderPre num2str(cluster+1) SubFolder 'IdxMatch' num2str(embryo(e)) '.mat']);
load([FolderPre num2str(cluster+1) SubFolder 'IdxMatchCentroid' num2str(embryo(e)) '.mat']);

    if anisotropy==1
        color = newAnisotropy(IdxMatchCentroid); 
    else
        color = areas(IdxMatchCentroid); 
    end

node2 = node(IdxMatch,:);
conn = meshconn(faceRef(:,1:3),size(node2,1));
n1{1} = node2;
niter = 2;

    for iter = 2:niter+1
        n1{iter} = smoothsurf(n1{iter-1},[],conn,1,0.9,'laplacian');
    end
    
node2 = n1{end};
plotmesh(node,face)
hold on
plotmesh(nodeRef,faceRef)

myo.vertices = node2;
myo.faces = faceRef;

save([FolderPre num2str(cluster) '\Deformation\DynamicAtlas\' singlecolor num2str(embryo(e)) '.mat'],'color');

write_ply(node2,faceRef,[FolderPre num2str(cluster) SubFolder 'new' num2str(embryo(e)) '.ply']);
delete ([FolderPre num2str(cluster) SubFolder singlecolor num2str(embryo(e)) '.mtl']);
delete ([FolderPre num2str(cluster) SubFolder singlecolor num2str(embryo(e)) '.obj']);
obj_write_color(myo, [FolderPre num2str(cluster) SubFolder singlecolor num2str(embryo(e))],color);
end



clear;close all; clc
embryo = [31,5];
cluster = 3;
anisotropy = 0;
FolderPre = '\Code_data\2. IntegratingMultipleLiveImagesIntoAConsensusTemporalReference\Cluster';
SubFolder = '\Deformation\DynamicAtlas\';


for e = 1:size(embryo,2)
    if anisotropy==1
        name = 'colormeanAnisotropy';
        singlecolor = 'colorAni';
        defor = 'Anisotropysmooth';
    else
        name = 'colormean';
        singlecolor = 'colorGro';
        defor ='JCCsmooth';
    end
    
[node,face] = read_ply([FolderPre num2str(cluster) SubFolder 'registered' num2str(embryo(e)) '.ply']);
[nodeRef, faceRef] = read_ply([FolderPre num2str(cluster) SubFolder 'new22.ply']);
IdxMatch = knnsearch(node, nodeRef);
save([FolderPre num2str(cluster) SubFolder 'IdxMatch' num2str(embryos(e)) '.mat'], 'IdxMatch');

centroid = meshcentroid(node, face);
centroidRef = meshcentroid(nodeRef, faceRef);
IdxMatchCentroid = knnsearch(centroid, centroidRef);
save([FolderPre num2str(cluster) SubFolder 'IdxMatchCentroid' num2str(embryos(e)) '.mat'], 'IdxMatchCentroid');


[nodeRef,faceRef] = read_ply([FolderPre num2str(cluster) '\3D\Mapped\27.ply']);
load ([FolderPre num2str(cluster) '\Deformation\' defor num2str(embryo(e)) '.mat'])

node2 = node(IdxMatch,:);

conn = meshconn(faceRef(:,1:3),size(node2,1));
n1{1} = node2;
niter = 2;

    for iter = 2:niter+1
        n1{iter} = smoothsurf(n1{iter-1},[],conn,1,0.9,'laplacian');
    end
    
node2 = n1{end};

write_ply(node2,faceRef,[FolderPre num2str(cluster) SubFolder 'new' num2str(embryo(e)) '.ply']);

    if anisotropy==1
        color = newAnisotropy(IdxMatchCentroid); 
    else
        color = areas(IdxMatchCentroid); 
    end

myo.vertices = node2;
myo.faces = faceRef;

delete ([FolderPre num2str(cluster) SubFolder singlecolor num2str(embryo(e)) '.mtl']);
delete ([FolderPre num2str(cluster) SubFolder singlecolor num2str(embryo(e)) '.obj']);
obj_write_color(myo, [FolderPre num2str(cluster) SubFolder singlecolor num2str(embryo(e))],color);

colormerge(:,e) = color;

    if e == size(embryo,2)
        [node,face]=read_ply([FolderPre num2str(cluster) SubFolder 'new22.ply']);
        myo.vertices = node;
        myo.faces = face;
        load ([FolderPre num2str(cluster) '\Deformation\' defor  '16.mat'])
        load ([FolderPre num2str(cluster+1) SubFolder 'IdxMatch16.mat'])
        load ([FolderPre num2str(cluster+1) SubFolder 'IdxMatchCentroid16.mat']);
        
            if anisotropy==1
                colormerge(:,e+1) = newAnisotropy(IdxMatchCentroid);
            else
                colormerge(:,e+1) = newJCC(IdxMatchCentroid);
            end
            
        load ([FolderPre num2str(cluster) '\Deformation\' defor  '22.mat'])
        load ([FolderPre num2str(cluster+1) SubFolder '\IdxMatch22.mat'])
        load ([FolderPre num2str(cluster+1) SubFolder '\IdxMatchCentroid22.mat']);
        
            if anisotropy==1
                colormerge(:,e+2) = newAnisotropy(IdxMatchCentroid);
            else
                colormerge(:,e+2) = newJCC(IdxMatchCentroid);
            end
            
        load ([FolderPre num2str(cluster) '\Deformation\' defor  '124.mat'])
        load ([FolderPre num2str(cluster+1) SubFolder '\IdxMatch124.mat'])
        load ([FolderPre num2str(cluster+1) SubFolder '\IdxMatchCentroid124.mat']);   
        
            if anisotropy==1
                colormerge(:,e+3) = newAnisotropy(IdxMatchCentroid);
            else
                colormerge(:,e+3) = areas(IdxMatchCentroid);
            end
            
        colormean = mean(colormerge,2);
        save([FolderPre num2str(cluster) SubFolder name '.mat'],'colormean');
        delete ([FolderPre num2str(cluster) SubFolder name '.mtl']);
        delete ([FolderPre num2str(cluster) SubFolder name '.obj']);
        obj_write_color(myo, [FolderPre num2str(cluster) SubFolder name ],colormean);
        write_ply(node,face,[FolderPre num2str(cluster) SubFolder 'mean.ply']);
    end
end


%% Propagation Gr5
clear;close all; clc
embryo= [35,7];
cluster = 5;
anisotropy = 0;
FolderPre = '\Code_data\2. IntegratingMultipleLiveImagesIntoAConsensusTemporalReference\';
SubFolder = '\Deformation\DynamicAtlas\';

for e = 1:size(embryo,2)
    if anisotropy==1
        name = 'colormeanAnisotropy';
        singlecolor = 'colorAni';
        defor = 'Anisotropysmooth';
    else
        name = 'colormean';
        singlecolor = 'colorGro';
        defor ='JCCsmooth';
        name = 'SurfArea';
    end

[node,face] = read_ply([FolderPre num2str(cluster) '\3D\Mapped\' num2str(embryo(e)) '.ply']);
[nodeRef,faceRef] = read_ply([FolderPre num2str(cluster) '\3D\Mapped\27.ply']);
load([FolderPre num2str(cluster) '\Deformation\' defor num2str(embryo(e)) '.mat'])
load([FolderPre' num2str(cluster-1) SubFolder 'IdxMatch' num2str(embryo(e)) '.mat'])
load([FolderPre num2str(cluster-1) SubFolder 'IdxMatchCentroid' num2str(embryo(e)) '.mat']);

node2 = node(IdxMatch,:);

conn = meshconn(faceRef(:,1:3),size(node2,1));
n1{1} = node2;
niter = 2;

    for iter = 2:niter+1
        n1{iter} = smoothsurf(n1{iter-1},[],conn,1,0.9,'laplacian');
    end

node2 = n1{end};

    if anisotropy==1
        color = newAnisotropy(IdxMatchCentroid); 
    else
        color = areas(IdxMatchCentroid); 
    end

myo.vertices = node2;
myo.faces = faceRef;

write_ply(node2,faceRef,[FolderPre num2str(cluster) SubFolder 'new' num2str(embryo(e)) '.ply']);
delete ([FolderPre num2str(cluster) SubFolder singlecolor num2str(embryo(e)) '.mtl']);
delete ([FolderPre num2str(cluster) SubFolder singlecolor num2str(embryo(e)) '.obj']);
obj_write_color(myo, [FolderPre num2str(cluster) SubFolder singlecolor num2str(embryo(e))],color);

colormerge(:,e) = color;

    if e == size(embryo,2)
        if anisotropy==1
            load ([FolderPre num2str(cluster) '\Deformation\Anisotropysmooth27.mat'])
            colormerge(:,e+1) = newAnisotropy;  
        else
            load ([FolderPre num2str(cluster) '\Deformation\SurfArea27.mat'])
            colormerge(:,e+1) = areas;   
        end
    colormean = mean(colormerge,2);
    save([FolderPre num2str(cluster) SubFolder name '.mat'],'colormean');
    delete ([FolderPre num2str(cluster) SubFolder name '.mtl']);
    delete ([FolderPre num2str(cluster) SubFolder name '.obj']);
    obj_write_color(myo, [FolderPre num2str(cluster) SubFolder name ],colormean);
    write_ply(node2,faceRef,[FolderPre num2str(cluster) SubFolder 'mean.ply']);
    end
end


%% Propagation Gr6
clear;close all; clc
embryo= 27;
cluster = 6;
FolderPre = '\Code_data\2. IntegratingMultipleLiveImagesIntoAConsensusTemporalReference\';
SubFolder = '\Deformation\DynamicAtlas\';

for e = 1:size(embryo,2)
[node,face] = read_ply([FolderPre num2str(cluster) '\3D\Mapped\' num2str(embryo(e)) '.ply']);
myo.vertices = node;
myo.faces = face;

load ([FolderPre num2str(cluster) '\Deformation\Anisotropysmooth' num2str(embryo(e)) '.mat']);
colormean = newAnisotropy;
delete([FolderPre num2str(cluster) SubFolder 'colormeanAnisotropy.mtl'])
delete([FolderPre num2str(cluster) SubFolder 'colormeanAnisotropy.obj'])
obj_write_color(myo, [FolderPre num2str(cluster) SubFolder 'colormeanAnisotropy'],colormean)
save([FolderPre num2str(cluster) SubFolder 'colormeanAnisotropy.mat'],'colormean');

load ([FolderPre num2str(cluster) '\Deformation\JCCsmooth' num2str(embryo(e)) '.mat'])
colormean = newJCC;
delete ([FolderPre num2str(cluster) SubFolder 'colormean.mtl'])
delete ([FolderPre num2str(cluster) SubFolder 'colormean.obj'])
obj_write_color(myo, [FolderPre num2str(cluster) SubFolder 'colormean'],colormean)
save([FolderPre num2str(cluster) SubFolder 'colormean.mat'],'colormean');

write_ply(node,face,[FolderPre num2str(cluster) SubFolder 'mean.ply']);
end

clear;close all; clc
embryo= 24;
cluster = 6;
FolderPre = '\Code_data\2. IntegratingMultipleLiveImagesIntoAConsensusTemporalReference\';
SubFolder = '\Deformation\DynamicAtlas\';

for e = 1:size(embryo,2)
[node,face] = read_ply([FolderPre num2str(cluster) SubFolder 'registered' num2str(embryo(e)) '.ply']);
conn = meshconn(face(:,1:3),size(node,1));
n1{1} = node;
niter = 2;

    for iter = 2:niter+1
        n1{iter} = smoothsurf(n1{iter-1},[],conn,1,0.9,'laplacian');
    end
node = n1{end};

[nodeRef,faceRef] = read_ply([FolderPre num2str(cluster) '\3D\Mapped\27.ply']);
IdxMatch = knnsearch(node,nodeRef);
save([FolderPre num2str(cluster) SubFolder 'IdxMatch' num2str(embryo(e)) '.mat'],'IdxMatch');
centroid = meshcentroid(node,face);
centroidRef = meshcentroid(nodeRef,faceRef);
IdxMatchCentroid = knnsearch(centroid,centroidRef);
save([FolderPre num2str(cluster) SubFolder 'IdxMatchCentroid' num2str(embryo(e)) '.mat'],'IdxMatchCentroid');
end


%% Propagation Gr7
clear;close all; clc
embryo= 24;
cluster = 7;
anisotropy = 0;
FolderPre = '\Code_data\2. IntegratingMultipleLiveImagesIntoAConsensusTemporalReference\Cluster';
SubFolder = '\Deformation\DynamicAtlas\';

for e = 1:size(embryo,2)
    if anisotropy==1
        name = 'colormeanAnisotropy';
        singlecolor = 'colorAni';
        defor = 'Anisotropysmooth';
    else
    name = 'colormean';
    singlecolor = 'colorGro';
    defor ='JCCsmooth';
     end
    
[node,face] = read_ply([FolderPre num2str(cluster) '\3D\Mapped\' num2str(embryo(e)) '.ply']);
[nodeRef,faceRef]=read_ply([FolderPre num2str(cluster-1) '\3D\Mapped\27.ply']);

load([FolderPre' num2str(cluster) '\Deformation\' defor num2str(embryo(e)) '.mat'])
load([FolderPre num2str(cluster-1) SubFolder 'IdxMatch' num2str(embryo(e)) '.mat'])
load([FolderPre num2str(cluster-1) SubFolder 'IdxMatchCentroid' num2str(embryo(e)) '.mat']);

node2 = node(IdxMatch,:);

conn = meshconn(faceRef(:,1:3),size(node2,1));
n1{1} = node2;
niter = 2;

    for iter = 2:niter+1
        n1{iter} = smoothsurf(n1{iter-1},[],conn,1,0.9,'laplacian');
    end

node2 = n1{end};

write_ply(node2,faceRef,[FolderPre num2str(cluster) SubFolder 'new' num2str(embryo(e)) '.ply']);

    if anisotropy==1
        color = newAnisotropy(IdxMatchCentroid); 
    else
        color = areas(IdxMatchCentroid); 
    end

myo.vertices = node2;
myo.faces = faceRef;

delete ([FolderPre num2str(cluster) SubFolder singlecolor num2str(embryo(e)) '.mtl']);
delete ([FolderPre num2str(cluster) SubFolder singlecolor num2str(embryo(e)) '.obj']);
obj_write_color(myo, [FolderPre num2str(cluster) SubFolder singlecolor num2str(embryo(e))],color);

colormerge(:,e) = color;

    if e == size(embryo,2)
        [node,face] = read_ply([FolderPre num2str(cluster) SubFolder 'new24.ply']);
        myo.vertices = node;
        myo.faces = face;
        colormean = mean(colormerge,2);
        save([FolderPre num2str(cluster) SubFolder name '.mat'],'colormean');
        delete([FolderPre num2str(cluster) SubFolder name '.mtl'])
        delete([FolderPre num2str(cluster) SubFolder name '.obj'])
        obj_write_color(myo, [FolderPre num2str(cluster) SubFolder name ],colormean)
        write_ply(node,face,[FolderPre num2str(cluster) SubFolder 'mean.ply']);
      
    end
end


%% Propagation Gr8
clear;close all; clc
embryo= 24;
cluster = 8;
anisotropy = 0;
FolderPre = '\Code_data\2. IntegratingMultipleLiveImagesIntoAConsensusTemporalReference\Cluster';
SubFolder = '\Deformation\DynamicAtlas\';

for e = 1:size(embryo,2)
    if anisotropy==1
        name = 'colormeanAnisotropy';
        singlecolor = 'colorAni';
        defor = 'Anisotropysmooth';
    else
    name = 'colormean';
    singlecolor = 'colorGro';
    defor ='JCCsmooth';
    end
    
[node,face] = read_ply([FolderPre num2str(cluster) '\3D\Mapped\' num2str(embryo(e)) '.ply']);
[nodeRef,faceRef] = read_ply([FolderPre num2str(cluster-2) '\3D\Mapped\27.ply']);

load ([FolderPre num2str(cluster) '\Deformation\' defor num2str(embryo(e)) '.mat'])
load ([FolderPre num2str(cluster-2) SubFolder 'IdxMatch' num2str(embryo(e)) '.mat'])
load ([FolderPre num2str(cluster-2) SubFolder 'IdxMatchCentroid' num2str(embryo(e)) '.mat']);

node2 = node(IdxMatch,:);

conn = meshconn(faceRef(:,1:3),size(node2,1));
n1{1} = node2;
niter = 2;

    for iter = 2:niter+1
       n1{iter} = smoothsurf(n1{iter-1},[],conn,1,0.9,'laplacian');
    end

node2 = n1{end};
write_ply(node2,faceRef,[FolderPre num2str(cluster) SubFolder 'new' num2str(embryo(e)) '.ply']);

    if anisotropy==1
        color = newAnisotropy(IdxMatchCentroid); 
    else
        color = areas(IdxMatchCentroid); 
    end

myo.vertices = node2;
myo.faces = faceRef;

delete([FolderPre num2str(cluster) SubFolder singlecolor num2str(embryo(e)) '.mtl']);
delete([FolderPre num2str(cluster) SubFolder singlecolor  num2str(embryo(e)) '.obj']);
obj_write_color(myo, [FolderPre num2str(cluster) SubFolder singlecolor  num2str(embryo(e))],color);

colormerge(:,e) = color;

    if e == size(embryo,2)
        [node,face] = read_ply([FolderPre num2str(cluster) SubFolder 'new24.ply']);
        myo.vertices = node;
        myo.faces = face;
        colormean = mean(colormerge,2);
        save([FolderPre num2str(cluster) SubFolder name '.mat'],'colormean');
        delete ([FolderPre num2str(cluster) SubFolder name '.mtl']);
        delete ([FolderPre num2str(cluster) SubFolder name '.obj']);
        obj_write_color(myo,[FolderPre num2str(cluster) SubFolder name ],colormean);
        write_ply(node,face,[FolderPre num2str(cluster) SubFolder '\mean.ply']);
          end
end


clear;close all; clc
embryo = 12;
cluster = 8;
FolderPre = '\\2. IntegratingMultipleLiveImagesIntoAConsensusTemporalReference\Cluster';
SubFolder = '\Deformation\DynamicAtlas\';

for e = 1:size(embryo,2)
[node,face] = read_ply([FolderPre num2str(cluster) SubFolder 'registered' num2str(embryo(e)) '.ply']);
conn = meshconn(face(:,1:3),size(node,1));
n1{1} =node;
niter = 2;

    for iter = 2:niter+1
        n1{iter}=smoothsurf(n1{iter-1},[],conn,1,0.9,'laplacian');
    end
    
node = n1{end};

[nodeRef,faceRef] = read_ply([FolderPre num2str(cluster) SubFolder 'new24.ply']);
IdxMatch = knnsearch(node,nodeRef);

save([FolderPre num2str(cluster) SubFolder 'IdxMatch' num2str(embryo(e)) '.mat'],'IdxMatch');
centroid = meshcentroid(node,face);
centroidRef = meshcentroid(nodeRef,faceRef);
IdxMatchCentroid = knnsearch(centroid,centroidRef);
save([FolderPre num2str(cluster) SubFolder 'IdxMatchCentroid' num2str(embryo(e)) '.mat'],'IdxMatchCentroid');
end

clear;close all; clc
embryo = 12;
cluster = 8;
anisotropy = 0;
FolderPre = '\Code_data\2. IntegratingMultipleLiveImagesIntoAConsensusTemporalReference\Cluster';
SubFolder = '\Deformation\DynamicAtlas\';


for e = 1:size(embryo,2)
    if anisotropy==1
        name = 'colormeanAnisotropy';
        singlecolor = 'colorAni';
        defor = 'Anisotropysmooth';
    else
    name = 'colormean';
    singlecolor = 'colorGro';
    defor ='JCCsmooth';
    end
    
[node,face] = read_ply([FolderPre num2str(cluster) SubFolder 'registered' num2str(embryo(e)) '.ply']);
[nodeRef,faceRef] = read_ply([FolderPre num2str(cluster-3) '\3D\Mapped\27.ply']);
load ([FolderPre num2str(cluster) '\Deformation\' defor num2str(embryo(e)) '.mat'])
load ([FolderPre num2str(cluster) SubFolder 'IdxMatch' num2str(embryo(e)) '.mat'])
load ([FolderPre num2str(cluster) SubFolder 'IdxMatchCentroid' num2str(embryo(e)) '.mat']);

node2 = node(IdxMatch,:);

conn = meshconn(faceRef(:,1:3),size(node2,1));
n1{1} = node2;
niter = 2;

    for iter = 2:niter+1
        n1{iter} = smoothsurf(n1{iter-1},[],conn,1,0.9,'laplacian');
    end
    
node2 = n1{end};

write_ply(node2,faceRef,[FolderPre num2str(cluster) SubFolder 'new' num2str(embryo(e)) '.ply']);

    if anisotropy==1
        color = newAnisotropy(IdxMatchCentroid); 
    else
        color = areas(IdxMatchCentroid); 
    end

myo.vertices = node2;
myo.faces = faceRef;

delete ([FolderPre num2str(cluster) SubFolder singlecolor num2str(embryo(e)) '.mtl']);
delete ([FolderPre num2str(cluster) SubFolder singlecolor num2str(embryo(e)) '.obj']);
obj_write_color(myo, [FolderPre num2str(cluster) SubFolder singlecolor num2str(embryo(e))],color);
colormerge(:,e) = color;

    if e == size(embryo,2)
        [node,face] = read_ply([FolderPre num2str(cluster) SubFolder 'new24.ply']);
        myo.vertices = node;
        myo.faces = face;
        load ([FolderPre num2str(cluster) '\Deformation\' defor  '24.mat'])
        load ([FolderPre num2str(cluster-2) SubFolder 'IdxMatch24.mat'])
        load ([FolderPre num2str(cluster-2) SubFolder 'IdxMatchCentroid24.mat']);
        
            if anisotropy==1
                colormerge(:,e+1) = newAnisotropy(IdxMatchCentroid);
            else
                colormerge(:,e+1) = neJCC(IdxMatchCentroid);
            end
            
        colormean = mean(colormerge,2);
        save([FolderPre num2str(cluster) SubFolder name '.mat'],'colormean');
        delete ([FolderPre num2str(cluster) SubFolder name '.mtl']);
        delete ([FolderPre num2str(cluster) SubFolder name '.obj']);
        obj_write_color(myo, [FolderPre num2str(cluster) SubFolder name ],colormean);
        write_ply(node,face,[FolderPre num2str(cluster) SubFolder 'mean.ply']);
    end
end

%% 3. Sum deformation from Gr2
clear;close all; clc

FolderPre = '\Code_data\2. IntegratingMultipleLiveImagesIntoAConsensusTemporalReference\Cluster';
SubFolder = '\Deformation\DynamicAtlas\';

Group= [3,4,5,6,7,8,9];
[nodeRef,faceRef] = read_ply([FolderPre num2str(Group(4)) '\3D\Mapped\27.ply']);

colorsum = ones(size(faceRef,1),1); % number of faces in embryo ref 
anisotropy = 0;

for e = 1:size(Group,2)
    if anisotropy==1
        name = 'colormeanAnisotropy';
        sum = 'SumAnisotropy';
    else
        name = 'colormean';
        sum = 'SumJCC';
    end
    
load ([FolderPre num2str(Group(e)) SubFolder name '.mat']);

colorsum = colormean.*colorsum;

[node,face]=read_ply([FolderPre num2str(Group(e)) SubFolder 'mean.ply']);
 
myo.vertices = node;
myo.faces= face ;

save([FolderPre num2str(Group(e)) SubFolder sum '.mat'],'colorsum');
delete ([FolderPre num2str(Group(e)) SubFolder sum '.mtl']);
delete ([FolderPre num2str(Group(e)) SubFolder sum '.obj']);
obj_write_color(myo,[FolderPre num2str(Group(e)) SubFolder sum ''],colorsum);
end

%% 4. Sum deformation into ATLAS from Gr2
clear; close all; clc;

% Define class IDs and corresponding indices
group = [3, 4, 5, 6, 7, 8, 9];
index = [22, 27, 27, 27, 24, 24, 12];

% Base folder path for data
FolderPre = '\Code_data\2. IntegratingMultipleLiveImagesIntoAConsensusTemporalReference\Cluster';
SubFolder = '\Deformation\DynamicAtlas\';


for e = 1:size(group,2)
        % Set output path
    outpath = ['V:\CellReportsMethods_2024\Data\Figure4_Quantifying_Tissue_Deformation_During_Early_Morphogenesis\Gr' num2str(group(e)) '\Cumulative\'];

    % Prepare and create output folder if it doesn't exist
    if ~exist(outpath, 'dir')
        mkdir(outpath);
    end

    % Loop through the two types: JCC (Growth) and Anisotropy
    for typeIdx = 1:2
        if typeIdx == 1
            name = 'ATLASsumJCC';
            sumFile = 'SumJCC';
            deformationType = 'Growth';
        else
            name = 'ATLASsumAnisotropy';
            sumFile = 'SumAnisotropy';
            deformationType = 'Anisotropy';
        end

    load([FolderPre num2str(group(e)) '\Deformation\DynamicAtlas\' sumFile '.mat']);


    if index(e)==27
        load([FolderPre num2str(group(e)) '\VertexCorrespondence\Embryo' num2str(index(e)) '\IdxCUT.mat']);
        load([FolderPre num2str(group(e)) '\VertexCorrespondence\Embryo' num2str(index(e)) '\IdxMatch.mat']);
    else
        [nodeA,faceA] = read_ply([FolderPre num2str(group(e)) '\3D\Atlas\' num2str(index(e)) '.ply']);
        [node,face] = read_ply([FolderPre num2str(group(e)) '\Deformation\DynamicAtlas\mean.ply']);

        % [tform, movingReg, rmse] = pcregistericp(pointCloud(node), pointCloud(nodeA));
        % node = movingReg.Location;
        % IdxMatch = knnsearch(node,nodeA);
        % 
        % write_ply(node(IdxMatch,:),faceA,[FolderPre num2str(classe(e)) '\Deformation\DynamicAtlas\meanRegistered.ply']);
        % this for the cuttin part


        [nodeC,faceC] = read_ply([FolderPre num2str(group(e)) '\Deformation\DynamicAtlas\Atlas_Cut.ply']);
        centroid = meshcentroid(nodeC,faceC);
        centroid1 = meshcentroid(nodeA,faceA);
        [IdxCUT,~] = knnsearch(centroid1,centroid);   
        % find the metching part
        [tform, movingReg, rmse] = pcregistericp(pointCloud(node), pointCloud(nodeC));
        node = movingReg.Location;
        centroid = meshcentroid(node,face);
        centroid1 = meshcentroid(nodeA,faceA);
        [IdxMatch,D] = knnsearch(centroid,centroid1(IdxCUT,:));   
    end

[node,face] = read_ply(['\\Source_Data\Figure 3\3A\Gr' num2str(group(e)) '.ply']);
  % Initialize color data
        color1 = zeros(size(face, 1), 1);
        color1(IdxCUT, 1) = colorsum(IdxMatch, 1);
        
        % Set color to zero where no match was found
        ind = any(color1 == 0, 2);
        color1(ind, :) = 0;
        
        % Generate RGB data for visualization
        rgbData = dataToRGB(color1, 1, typeIdx, 'mean');
        
        % Save the mesh with color information
        writeMesh_plyModify([outpath name '.ply'], node, face, color1, rgbData);
    end
clearvars -except group index anisotropy mincol maxcol maps FolderPre
end

