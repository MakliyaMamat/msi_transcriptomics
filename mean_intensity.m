% Regions as 11x2 cell array and 'energy' as 213x25469 matrix
numRegions = size(t_values, 1); 
numGenes = size(energy, 2);    

% Extract the structure names from structInfo (assuming names are in the 4th column)
structure_names = string(structInfo{:, 4});  % Convert to string array for easy comparison
structure_names = strtrim(lower(structure_names));  % Trim and convert to lowercase

% Initialize a matrix to store the mean expression intensity for each region
mean_gene_energy = NaN(numRegions, numGenes);

% Loop over each region
for i = 1:numRegions
    % Get the structure names for the current region
    region_structures = Regions{i, 2};
    
    % Convert structure names to lowercase for comparison
    region_structures = lower(strtrim(region_structures));
    
    % Find the indices of structures that match the current region
    structure_indices = [];
    for j = 1:length(region_structures)
        % Use 'contains' to match structures robustly
        idx = find(contains(structure_names, region_structures(j), 'IgnoreCase', true));
        structure_indices = [structure_indices; idx];
    end
    
    structure_indices = unique(structure_indices);   % Remove duplicates
    structure_indices = structure_indices(~isnan(structure_indices) & structure_indices > 0); % Remove NaN
    
    % Check if there are valid indices for region
    if ~isempty(structure_indices)
        % Extract the energy values for the structures of the current region
        region_energy_values = energy(structure_indices, :);

        % Compute the mean energy across the structures
        mean_gene_energy(i, :) = nanmean(region_energy_values, 1);
    else
        warning(['No valid indices found for region: ', Regions{i, 1}]);
    end
end

% 'mean_gene_energy' saved as 11, 25469
