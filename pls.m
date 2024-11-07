% Load gene expression data (11x25469 matrix)
filename = 'ms_cleaned_data_TET2.xlsx';  
ms_cleaned_data = readmatrix(filename);  
filename = 'Tet2pvalues.xlsx';  
p_values = readmatrix(filename);  
filename = 'TET2_tvalues.xlsx';  
t_values = readmatrix(filename);  

p_threshold = 0.05;
significant_metabolites = any(p_values < p_threshold, 1);  % Identify significant metabolites
significant_metabolite_indices = find(significant_metabolites);  % Indices of significant metabolites
t_values_significant = t_values(:, significant_metabolites);
p_values_significant = p_values(:, significant_metabolites);  

min_p_values = min(p_values_significant, [], 1);  % Take the minimum across the rows (regions)
[~, most_significant_index] = min(min_p_values);
p_values_most_significant = p_values_significant(:, most_significant_index);  % Corresponding p-values
disp(['The most significant metabolite is at position ', num2str(most_significant_index)]);
disp(['Minimum p-value for this metabolite is: ', num2str(min_p_values(most_significant_index))]);

% Select the metabolite of interest by its index in the original set
original_metabolite_index = 41;  
selected_metabolite_position = find(significant_metabolite_indices == original_metabolite_index);
t_values_col = t_values_significant(:, selected_metabolite_position);  % Extract the column for the selected metabolite
t_values_col = table2array(t_values_col);  % Convert

% Filter out NaNs from t-values and corresponding rows from gene expression data
valid_indices_tvalues = ~isnan(t_values_col);  % Find indices without NaNs
Tvalues_filtered = t_values_col(valid_indices_tvalues);  % Filter t-values based on valid indices
gene_filtered = aggregated_gene_energy(valid_indices_tvalues, :);  % Filter gene expression data based on valid tvalue indices

% Remove columns (genes) in gene_filtered that contain any NaN
valid_indices_gene = all(~isnan(gene_filtered), 1);  % Find columns without NaNs in gene expression data
gene_filtered = gene_filtered(:, valid_indices_gene);  % Keep only valid columns (genes) in the filtered gene expression data

% Filter the gene names the same way as X_filtered
gene_names = geneInfo{:, 1};  % Assuming you have the original gene names in geneInfo
filtered_gene_names = gene_names(valid_indices_gene);  % Keep only the gene names that correspond to the valid columns

X = gene_filtered
Y = Tvalues_filtered; % Extract t-values (as one column vector)

% z-score normalization of gene expression values and t-values
X = zscore(gene_filtered);
Y = zscore(Y);

% Perform full PLS and plot variance in Y explained by top 8 components
[XL, YL, XS, YS, BETA, PCTVAR, MSE, stats] = plsregress(X, Y);
dim = 8;
plot(1:dim, cumsum(100 * PCTVAR(2, 1:dim)), '-o', 'LineWidth', 1.5, 'Color', [140/255, 0, 0]);
set(gca, 'FontSize', 14);
xlabel('Number of PLS components', 'FontSize', 14);
ylabel('Percent Variance Explained in Y', 'FontSize', 14);
grid on;
saveas(gca, 'PLS_ExplainedVariance.png');

dim=2;
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim); % no need to do this but it keeps outputs tidy

%%% plot correlation of PLS component 1 with t-statistic (from Cobre as an example):
figure
plot(XS(:,1), Y, 'r.')  
xlabel('PLS Component 1 Scores', 'FontSize', 10);  
ylabel('Standardized T values', 'FontSize', 10);  
title('Correlation of PLS Component 1 with T-value', 'FontSize', 10);
grid on
saveas(gca, 'Corr_PLS1_Tvalue.png');
figure
plot(XS(:,2), Y, 'r.')  
xlabel('PLS Component 2 Scores', 'FontSize', 10);  
ylabel('Standardized T values', 'FontSize', 10);  
title('Correlation of PLS Component 2 with T-value', 'FontSize', 10);
grid on
saveas(gca, 'Corr_PLS2_Tvalue.png');

% permutation testing to assess significance of PLS result as a function of
% the number of components (dim) included:
rep=1000;
for dim=1:8
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim);
temp=cumsum(100*PCTVAR(2,1:dim));
Rsquared = temp(dim);
    for j=1:rep
        %j
        order=randperm(size(Y,1));
        Yp=Y(order,:);

        [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Yp,dim);

        temp=cumsum(100*PCTVAR(2,1:dim));
        Rsq(j) = temp(dim);
    end
dim
R(dim)=Rsquared
p(dim)=length(find(Rsq>=Rsquared))/rep
end

figure
plot(1:dim, p,'ok','MarkerSize',8,'MarkerFaceColor','b');
xlabel('Number of PLS components','FontSize',10);
ylabel('p-value','FontSize',10);
grid on
filename = 'PLS_components_pvalue.png';  % Define the filename
saveas(gcf, filename);  % Save the current figure as a PNG

dim=2;
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim);

% PLS in 2 dimensions 
dim=2;
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim);

%store regions' IDs and weights in descending order of weight for both components:
[R1,p1]=corr([XS(:,1),XS(:,2)],Y);

%align PLS components with desired direction for interpretability 
if R1(1,1)<0  
    stats.W(:,1)=-1*stats.W(:,1);
    XS(:,1)=-1*XS(:,1);
end
if R1(2,1)<0 
    stats.W(:,2)=-1*stats.W(:,2);
    XS(:,2)=-1*XS(:,2);
end

geneindex = 1:length(filtered_gene_names);
[PLS1w,x1] = sort(stats.W(:,1),'descend');
PLS1ids=filtered_gene_names(x1);
geneindex1=geneindex(x1);
[PLS2w,x2] = sort(stats.W(:,2),'descend');
PLS2ids=filtered_gene_names(x2);
geneindex2=geneindex(x2);

%define variables for storing the (ordered) weights from all bootstrap runs
PLS1weights=[];
PLS2weights=[];

bootnum=1000; %number of bootstrap iterations:

%start bootstrap
for i=1:bootnum
    i
    myresample = randsample(size(X,1),size(X,1),1);
    res(i,:)=myresample; %store resampling out of interest
    Xr=X(myresample,:); 
    Yr=Y(myresample,:); 
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(Xr,Yr,dim); % perform PLS for resampled data
      
    temp=stats.W(:,1);% extract PLS1 weights
    newW=temp(x1); % order the newly obtained weights the same way as initial PLS 
    if corr(PLS1w,newW)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
        newW=-1*newW;
    end
    PLS1weights=[PLS1weights,newW]; % store (ordered) weights from this bootstrap run
    
    temp=stats.W(:,2);
    newW=temp(x2); 
    if corr(PLS2w,newW)<0 
        newW=-1*newW;
    end
    PLS2weights=[PLS2weights,newW]; 
end

% get standard deviation of weights from bootstrap runs
PLS1sw=std(PLS1weights');
PLS2sw=std(PLS2weights');

% get bootstrap weights
temp1=PLS1w./PLS1sw';
temp2=PLS2w./PLS2sw';

% order bootstrap weights (Z) and names of regions
[Z1 ind1]=sort(temp1,'descend');
PLS1=PLS1ids(ind1);
geneindex1=geneindex1(ind1);
[Z2 ind2]=sort(temp2,'descend');
PLS2=PLS2ids(ind2);
geneindex2=geneindex2(ind2);

% save
% bootstrapped ordered list of genes
fid1 = fopen('PLS1_geneWeights.csv','w')
for i=1:length(filtered_gene_names)
  fprintf(fid1,'%s, %d, %f\n', PLS1{i}, geneindex1(i), Z1(i));
end
fclose(fid1)

fid2 = fopen('PLS2_geneWeights.csv','w')
for i=1:length(filtered_gene_names)
  fprintf(fid2,'%s, %d, %f\n', PLS2{i},geneindex2(i), Z2(i));
end
fclose(fid2)

filename_pvalues = 'PLS_PermutationTesting_Pvalues.xlsx';
p_table = array2table(p', 'VariableNames', {'P_value'});  % Convert to table and transpose p-values
writetable(p_table, filename_pvalues);

% Save PLS1 and PLS2 t-values as Excel files
filename_pls1 = 'PLS1_Tvalues.xlsx';
filename_pls2 = 'PLS2_Tvalues.xlsx';
PLS1_t_table = array2table(XS(:, 1), 'VariableNames', {'PLS1_Tvalues'});  % PLS1 t-values
PLS2_t_table = array2table(XS(:, 2), 'VariableNames', {'PLS2_Tvalues'});  % PLS2 t-values
writetable(PLS1_t_table, filename_pls1);
writetable(PLS2_t_table, filename_pls2);
disp('PLS1 and PLS2 t-values saved successfully.');

% Save PLS1 and PLS2 ROI scores
xlswrite('PLS1_ROIscores.csv',XS(:,1));
xlswrite('PLS2_ROIscores.csv',XS(:,2));
