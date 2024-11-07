% Extract the 'energy' field from GeneExpData (213x19419 matrix)
energy = GeneExpData.gene_energy;

% Extract the structure names from structInfo (names are in the 4th column)
structure_names = string(structInfo{:, 4});  
structure_names = strtrim(lower(structure_names));  % Remove any leading/trailing whitespace and convert to lowercase

% Define regions and their corresponding structures
Regions = {
    'Anterior Olfactory', ["Accessory olfactory bulb", "Anterior olfactory nucleus", "Main olfactory bulb"];
    'Basal Forebrain', ["Nucleus accumbens", "Substantia innominata", "Bed nuclei of the stria terminalis", "Magnocellular nucleus"];
    'Caudate Putamen', ["Caudoputamen", "Fundus of striatum"];
    'Cerebellum', ["Central lobule", "Culmen", "Dentate nucleus", "Flocculus", "Pyramus (VIII)", "Nodulus (X)", "Simple lobule", "Ansiform lobule", "Interposed nucleus", "Paramedian lobule", "Paraflocculus"];
    'Cortex', ["Posterior parietal association areas", "Anterior cingulate area, dorsal part", "Primary motor area", "Secondary motor area", "Primary somatosensory area, barrel field", "Primary visual area", "Retrosplenial area, dorsal part"];
    'Hippocampus', ["Field CA1", "Field CA2", "Field CA3", "Subiculum, dorsal part", "Subiculum, ventral part", "Dentate gyrus", "Parasubiculum", "Entorhinal area, lateral part", "Entorhinal area, medial part, dorsal zone", "Postsubiculum", "Presubiculum"];
    'Midbrain', ["Red nucleus", "Substantia nigra, compact part", "Superior colliculus, motor related", "Superior colliculus, sensory related", "Midbrain reticular nucleus", "Ventral tegmental area"];
    'Olfactory Bulbs', ["Accessory olfactory bulb", "Main olfactory bulb", "Nucleus of the lateral olfactory tract"];
    'Pons & Medulla', ["Pontine reticular nucleus", "Superior olivary complex", "Parabrachial nucleus", "Medullary reticular nucleus, dorsal part", "Medullary reticular nucleus, ventral part"];
    'Thalamus', ["Lateral dorsal nucleus of thalamus", "Mediodorsal nucleus of thalamus", "Anteroventral nucleus of thalamus", "Reticular nucleus of the thalamus", "Ventral anterior-lateral complex of the thalamus"];
    'Ventral Striatum', ["Nucleus accumbens", "Olfactory tubercle"];
};

% Initialize structure array to store filtered indices
filtered_indices = struct('Region', {}, 'Indices', {});

% Loop through each region to find and store indices of relevant structures
for i = 1:size(Regions, 1)
    region_name = Regions{i, 1};
    region_structures = strtrim(lower(Regions{i, 2})); % Convert structures to lowercase
    region_indices = [];

    % Find indices of structures matching the names in region_structures
    for j = 1:length(region_structures)
        idx = find(strcmp(structure_names, region_structures(j)));
        region_indices = [region_indices; idx];  % Append indices
    end

    % Store the found indices in the struct
    filtered_indices(i).Region = region_name;
    filtered_indices(i).Indices = region_indices;
end

% Extract unique regions from the struct array
unique_regions = unique({filtered_indices.Region});

% Initialize an array to store aggregated energy values for each unique region
aggregated_energy = NaN(length(unique_regions), size(energy, 2));

% Compute the mean energy for each unique region using the filtered indices
for i = 1:length(unique_regions)
    region_name = unique_regions{i};
    
    % Find the indices in the filtered structure corresponding to this region
    idx = find(strcmp({filtered_indices.Region}, region_name));
    
    % Collect all indices associated with this region
    region_specific_indices = [];
    for k = 1:length(idx)
        region_specific_indices = [region_specific_indices; filtered_indices(idx(k)).Indices];
    end
    
    % Check if valid indices are found and compute the mean
    if ~isempty(region_specific_indices)
        aggregated_energy(i, :) = nanmean(energy(region_specific_indices, :), 1);
    else
        warning(['No valid indices found for region: ', region_name]);
    end
end

% 'aggregated_energy' saved as 11x19419
