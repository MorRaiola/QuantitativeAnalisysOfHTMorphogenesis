% Function to save deformation data
function saveDeformationData(folder, embryoId, JCC, anisotropy, globalAutovet, barycenter2)
    save(fullfile(folder, ['JCC' num2str(embryoId) '.mat']), 'JCC');
    save(fullfile(folder, ['Anisotropy' num2str(embryoId) '.mat']), 'anisotropy');
    save(fullfile(folder, ['GlobalautovetCC' num2str(embryoId) '.mat']), 'globalAutovet');
    save(fullfile(folder, ['baricentro2' num2str(embryoId) '.mat']), 'barycenter2');
end