function registerGroup(embryo, reference, cluster,Folder,SubFolder)
Folder = [Folder cluster SubFolder];
    for e = 1:numel(embryo)
        % Load necessary data
        [node, face] = read_ply([Folder 'resize' num2str(embryo(e)) '.ply']);
        trasf = load([Folder num2str(embryo(e)) '.mat']).res;
        F = mirt3D_F(trasf.okno);
        imagRef = loadtiff([Folder reference '.tif']);
        imagTest = loadtiff([Folder 'resize' num2str(embryo(e)) '.tif']);

        % Node boundary correction
        node = max(min(node, [size(imagRef,2) size(imagRef,1) size(imagRef,3)]), 1);

        % Register nodes
        point1 = registerNodes(node, imagRef, imagTest, trasf);

        % Save registered nodes and faces
        write_ply(point1, face, [Folder 'registered' num2str(embryo(e)) '.ply']);
        [newnode2, newface2] = surfreorient(point1, face);
        write_ply(newnode2, newface2, [Folder 'registered' num2str(embryo(e)) '.ply']);

        % Smooth surface
        smoothAndSaveSurface(point1, face, Folder, embryo(e));

        clearvars -except embryo reference cluster
    end
end

function point1 = registerNodes(node, imagRef, imagTest, trasf)
    point1r = cell(size(node,1),1); point1c = cell(size(node,1),1); point1v = cell(size(node,1),1);

    tic
    for i = 1:size(node,1)
        cellIM = zeros(size(imagRef,1), size(imagRef,2));
        cellIM(int16(node(i,2)), int16(node(i,1))) = 1;
        [r,c] = find(cellIM);

        [r1,c1] = find(imagTest(:,:,int16(node(i,3))) ~= 0);
        [Idx, ~] = knnsearch([r1,c1], [r,c]);

        cells = zeros(size(imagTest));
        cells(int16(r1(Idx)), int16(c1(Idx)), int16(node(i,3))) = 1;
        rim = mirt3D_transform(double(cells), trasf);

        [r, c, v] = ind2sub(size(rim), find(rim ~= 0));
        [r1, c1, v1] = ind2sub(size(imagRef), find(imagRef ~= 0));
        [Idx, D] = knnsearch([r1,c1,v1], [r,c,v]);

        if ~isempty(D)
            pos = find(D == min(D), 1);
            point1r{i,1} = r1(Idx(pos));
            point1c{i,1} = c1(Idx(pos));
            point1v{i,1} = v1(Idx(pos));
        end
    end
    toc

    % Replacing zeros with nearest neighbors
    point1 = cell2mat([point1r point1c point1v]);
    [x, ~] = find(point1 == 0);
    [Idx, ~] = knnsearch(node, node(x,:));
    point1(x,:) = point1(Idx,:);
end

function smoothAndSaveSurface(node, face, Folder, embryo)
    conn = meshconn(face(:,1:3), size(node,1));
    n1{1} = node;
    niter = 4;

    for iter = 2:niter+1
        n1{iter} = smoothsurf(n1{iter-1}, [], conn, 1, 0.9, 'laplacian');
    end

    write_ply(n1{end}, face, [Folder 'registered' num2str(embryo) '.ply']);
end

