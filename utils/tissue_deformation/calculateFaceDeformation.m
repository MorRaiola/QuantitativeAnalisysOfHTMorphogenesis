function [F, JCC, anisotropy, globalAutovet] = calculateFaceDeformation(faceNodes1, faceNodes2, localSystem1, localSystem2, xv, yv, zv)
    % ** System 1 **: Compute translation and rotation matrix A1
    trasl1 = faceNodes1 - faceNodes1(1, :); % Translation to local coordinate system
    A1 = [xv; yv; zv] / [localSystem1(1,1:3);localSystem1(1,4:6);localSystem1(1,7:9)]; % Rotation matrix
    nodesNewCC1 = round(trasl1 * A1, 3);
    
    if all(nodesNewCC1(:,3) == 0)
        nodesNewCC1(:,3) = []; % Remove third column if it's all zeros
    end

    nodesNewCC1(1,nodesNewCC1(1,:) == 0) = 0.0001;
    nodesNewCC1(2,nodesNewCC1(2,:) == 0) = 0.0002;
    nodesNewCC1(3,nodesNewCC1(3,:) == 0) = 0.0003;

    % Compute edges e1 and e2 for system 1
    e1_1 = nodesNewCC1(2,:) - nodesNewCC1(1,:);
    e2_1 = nodesNewCC1(3,:) - nodesNewCC1(1,:);
    % Check for zero values and replace with 0.001
    R = [e1_1' e2_1']; % Matrix for system 1

    % ** System 2 **: Compute translation and rotation matrix A2
    trasl2 = faceNodes2 - faceNodes2(1, :); % Translation to local coordinate system
    A2 = [xv; yv; zv] / [localSystem2(1,1:3);localSystem2(1,4:6);localSystem2(1,7:9)]; % Rotation matrix
    nodesNewCC2 = round(trasl2 * A2, 3);
    
    if all(nodesNewCC2(:,3) == 0)
        nodesNewCC2(:,3) = []; % Remove third column if it's all zeros
    end

    nodesNewCC2(1,nodesNewCC2(1,:) == 0) = 0.001;
    nodesNewCC2(2,nodesNewCC2(2,:) == 0) = 0.002;
    nodesNewCC2(3,nodesNewCC2(3,:) == 0) = 0.003;
    % Compute edges e1 and e2 for system 2
    e1_2 = nodesNewCC2(2,:) - nodesNewCC2(1,:);
    e2_2 = nodesNewCC2(3,:) - nodesNewCC2(1,:);
    % Check for zero values and replace with 0.001
    T = [e1_2' e2_2']; % Matrix for system 2
    
    % Calculate the deformation gradient
    F = T / R;  % Deformation gradient

    % Calculate growth rate (JCC)
    JCC = det(F);

    % Compute Cauchy strain tensor
    CauchyStrain = F' * F;
    [evecs, evals]  = eig(CauchyStrain); % Eigenvalues of the strain tensor

    % Determine anisotropy
    anisotropy = sqrt(max(evals) / min(evals)); % Assuming evals(2) is the larger eigenvalue

    % Calculate the global autovectors
    % ordering the eigenvalues and eigenvectors in min max 
    [~,posmax] = max(diag(evals));
    [~,posmin] = min(diag(evals));
    eigvecs = [evecs(:,posmin)'; evecs(:,posmax)'];
    eigvecs(:,3) = 0;
    eigvecs(3,:) = 0;
    eigvecs(3,3) = 0;
    globalAutovet = eigvecs / A2;  % Transform eigenvectors to the global coordinate system
end
