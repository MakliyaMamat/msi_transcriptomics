% Load the MS data
[~, txt, raw] = xlsread('tet1_ms_data.xlsx');  

% Extract the sample names
sample_names = txt(2:end, 1);  

% Extract the intensity values
ms_intensity = cell2mat(raw(2:end, 2:end));  % Assuming this is numeric data

% Define the brain regions of interest
regions = {'AON', 'BFb', 'CPu', 'Cb', 'Ct', 'Hp', 'Mb', 'OB', 'Pn&Md', 'Th', 'VSa'};

% Initialize arrays to store the indices for TET and WD
TET_indices = cell(length(regions), 1);
WD_indices = cell(length(regions), 1);

% Find TET and WD sample indices
for r = 1:length(regions)
    region = regions{r};
    
    % Construct patterns 
    TET_pattern = ['AAV-cre-EGFP.*-' region '-'];
    WD_pattern = ['AAV-EGFP.*-' region '-'];

    TET_samples = find(~cellfun('isempty', regexp(sample_names, TET_pattern)));
    TET_indices{r} = TET_samples;
    WD_samples = find(~cellfun('isempty', regexp(sample_names, WD_pattern)));
    WD_indices{r} = WD_samples;
end

% Store the t-values (11 regions x 42 metabolites)
t_values = NaN(length(regions), size(ms_intensity, 2));

% Calculate t-values for the metabolites
for r = 1:length(regions)
    % indices for TET and WD
    tet_idx = TET_indices{r};
    wd_idx = WD_indices{r};
    
    % Check if there are enough samples to perform the t-test
    if ~isempty(tet_idx) && ~isempty(wd_idx)
        tet_intensities = ms_intensity(tet_idx, :);
        wd_intensities = ms_intensity(wd_idx, :);

        % Perform the t-test for each metabolite
        for m = 1:size(ms_intensity, 2)
            tet_values = tet_intensities(:, m);
            wd_values = wd_intensities(:, m);

            % Check if both groups have more than one valid (non-NaN) value
            if length(tet_values(~isnan(tet_values))) > 1 && length(wd_values(~isnan(wd_values))) > 1
                % Perform a two-sample t-test (use 'omitnan' to handle NaNs)
                [~, ~, ~, stats] = ttest2(tet_values, wd_values, 'Vartype', 'unequal');
                t_values(r, m) = stats.tstat;  % Store the t-value
            else
                warning(['Not enough valid data for metabolite ' num2str(m) ' in region: ' regions{r}]);
            end
        end
    else
        warning(['No valid TET or WD samples found for region: ' regions{r}]);
    end
end

