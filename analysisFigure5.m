%% Figure 5

%% 1.1 Correct CT for covariates: HCP

clc; clear; close all

load dataHCP
load data

clearvars -except dataHCP labelsGlasser*

dataHCP_Temp = dataHCP(:, {'ID' 'Age' 'IQ'});
CT_CorrectedHCP = dataHCP(:, labelsGlasserCT);
CT_CorrectedHCP = [dataHCP_Temp CT_CorrectedHCP];
CT_CorrectedHCP_Array = nan(height(CT_CorrectedHCP), width(CT_CorrectedHCP) - width(dataHCP_Temp));
clearvars dataHCP_Temp

measure = 'Thickness';
covariates = '~ Age + IQ';

measureOfInterest = strfind(CT_CorrectedHCP.Properties.VariableNames, measure);
measureOfInterestID = find(~cellfun(@isempty, measureOfInterest));

for parcelID = 1 : numel(measureOfInterestID)

    parcel = CT_CorrectedHCP.Properties.VariableNames{measureOfInterestID(parcelID)};

    fprintf('Correcting: %s\n', parcel)

    modelToBeFitted = sprintf('%s %s', parcel, covariates);
    model = fitlm(CT_CorrectedHCP, modelToBeFitted, 'RobustOpts', 'on');
    CT_CorrectedHCP_Array(:, parcelID) = model.Residuals.Raw;

end

CT_CorrectedHCP(:, labelsGlasserCT) = array2table(CT_CorrectedHCP_Array);

clearvars -except data* CT* labels*

save CT_CorrectedHCP

%% 1.2 Compute and plot CT G1: HCP

clc; clear; close all

load CT_CorrectedHCP

clearvars -except CT_CorrectedHCP*

% Compute and Fisher-transform CT covariance
CT_CovHCP = atanh(corr(CT_CorrectedHCP_Array));

% Nullify 'Inf' (along the diagonal)
CT_CovHCP(isinf(CT_CovHCP)) = 0;

% Create a gradient backbone using DME and normalised angle
gradientBone = GradientMaps('approach', 'dm', 'kernel', 'na');

% Fit the matrix to the backbone
gradientsCT_HCP = gradientBone.fit(CT_CovHCP);

% Derive the first 10 gradients
gradientsArrayCT_HCP = gradientsCT_HCP.gradients{1};

% Derive G1
G1_CT_HCP = gradientsArrayCT_HCP(:, 1);

save gradientsCT_HCP

% Plot G1
G1_CT_HCP_Conte69 = parcel_to_surface(G1_CT_HCP, 'glasser_360_conte69_fs_LR');
G1_CT_HCP_Conte69(G1_CT_HCP_Conte69 == 0) = (min(G1_CT_HCP) + max(G1_CT_HCP)) / 2;
figure('units', 'normalized', 'outerposition', [0.125 0 0.693 1])
title('CT G1: HCP', 'FontSize', 30)
axis off
plot_cortical(G1_CT_HCP_Conte69, 'surface_name', 'conte69', 'label_text', 'Loading');
colormap(flipud(brewermap([], 'PuOr')))
set(gca, 'FontSize', 25)

colourbars = findall(gcf, 'type', 'colorbar');
for cb = colourbars'
    set(cb, 'Ticks', [min(G1_CT_HCP) max(G1_CT_HCP)]);
    set(cb, 'TickLabels', {sprintf('%.3f', min(G1_CT_HCP)) sprintf('%.3f', max(G1_CT_HCP))});
end

%% 2.1 Correct CT for covariates: Total

clc; clear; close all

load data

clearvars -except data labelsGlasser*

dataTemp = data(:, {'URSI' 'PCL_Total' 'Status' 'Age' 'IQ'});
CT_Corrected = data(:, labelsGlasserCT);
CT_Corrected = [dataTemp CT_Corrected];
CT_CorrectedArray = nan(height(CT_Corrected), width(CT_Corrected) - width(dataTemp));
clearvars dataTemp

measure = 'Thickness';
covariates = '~ Age + IQ';

measureOfInterest = strfind(CT_Corrected.Properties.VariableNames, measure);
measureOfInterestURSI = find(~cellfun(@isempty, measureOfInterest));

for parcelID = 1 : numel(measureOfInterestURSI)

    parcel = CT_Corrected.Properties.VariableNames{measureOfInterestURSI(parcelID)};

    fprintf('Correcting: %s\n', parcel)

    modelToBeFitted = sprintf('%s %s', parcel, covariates);
    model = fitlm(CT_Corrected, modelToBeFitted, 'RobustOpts', 'on');
    CT_CorrectedArray(:, parcelID) = model.Residuals.Raw;

end

CT_Corrected(:, labelsGlasserCT) = array2table(CT_CorrectedArray);

clearvars -except data* CT* labels*

save CT_Corrected

%% 2.2 Compute and plot CT G1: Total

clc; clear; close all

load CT_Corrected

clearvars -except CT_Corrected*

% Compute and Fisher-transform CT covariance
CT_CovTotal = atanh(corr(CT_CorrectedArray));

% Nullify 'Inf' (along the diagonal)
CT_CovTotal(isinf(CT_CovTotal)) = 0;

% Create a gradient backbone using DME and normalised angle
gradientBone = GradientMaps('approach', 'dm', 'kernel', 'na');

% Fit the matrix to the backbone
gradientsCT_Total = gradientBone.fit(CT_CovTotal);

% Derive the first 10 gradients
gradientsArrayCT_Total = gradientsCT_Total.gradients{1};

% Derive G1
G1_CT_Total = gradientsArrayCT_Total(:, 1);

save gradientsCT_Total

% Plot G1
G1_CT_TotalConte69 = parcel_to_surface(G1_CT_Total, 'glasser_360_conte69_fs_LR');
G1_CT_TotalConte69(G1_CT_TotalConte69 == 0) = (min(G1_CT_Total) + max(G1_CT_Total)) / 2;
figure('units', 'normalized', 'outerposition', [0.125 0 0.693 1])
title('CT G1: Total', 'FontSize', 30)
axis off
plot_cortical(G1_CT_TotalConte69, 'surface_name', 'conte69', 'label_text', 'Loading');
colormap(flipud(brewermap([], 'PuOr')))
set(gca, 'FontSize', 25)

colourbars = findall(gcf, 'type', 'colorbar');
for cb = colourbars'
    set(cb, 'Ticks', [min(G1_CT_Total) max(G1_CT_Total)]);
    set(cb, 'TickLabels', {sprintf('%.3f', min(G1_CT_Total)) sprintf('%.3f', max(G1_CT_Total))});
end

%% 2.3 Compute and plot CT G1: Low

clc; clear; close all

load CT_Corrected
load gradientsCT_Total

clearvars -except CT_Corrected* gradients*

% Derive CT
CT_CorrectedArrayLow = CT_CorrectedArray(CT_Corrected.Status == 'Low', :);

% Compute and Fisher-transform CT covariance
CT_CovLow = atanh(corr(CT_CorrectedArrayLow));

% Nullify 'Inf' (along the diagonal)
CT_CovLow(isinf(CT_CovLow)) = 0;

% Create a gradient backbone using DME, normalised angle, and Procrustes 
% alignment
gradientBone = GradientMaps('approach', 'dm', 'kernel', 'na', 'alignment', 'Procrustes');

% Fit the matrix to the backbone with the total-sample gradients as the
% reference
gradientsCT_Low = gradientBone.fit(CT_CovLow, 'reference', gradientsArrayCT_Total);

% Derive the first 10 gradients
gradientsArrayCT_Low = gradientsCT_Low.aligned{1};

% Derive G1, aligned
G1_CT_Low = gradientsCT_Low.aligned{1}(:, 1);

save gradientsCT_Low

% Plot G1, aligned
G1_CT_LowConte69 = parcel_to_surface(G1_CT_Low, 'glasser_360_conte69_fs_LR');
G1_CT_LowConte69(G1_CT_LowConte69 == 0) = (min(G1_CT_Low) + max(G1_CT_Low)) / 2;
figure('units', 'normalized', 'outerposition', [0.125 0 0.693 1])
title('CT G1, aligned: Low', 'FontSize', 30)
axis off
plot_cortical(G1_CT_LowConte69, 'surface_name', 'conte69', 'label_text', 'Loading');
colormap(flipud(brewermap([], 'PuOr')))
set(gca, 'FontSize', 25)

colourbars = findall(gcf, 'type', 'colorbar');
for cb = colourbars'
    set(cb, 'Ticks', [min(G1_CT_Low) max(G1_CT_Low)]);
    set(cb, 'TickLabels', {sprintf('%.3f', min(G1_CT_Low)) sprintf('%.3f', max(G1_CT_Low))});
end

%% 2.4 Compute and plot CT G1: High

clc; clear; close all

load CT_Corrected
load gradientsCT_Total

clearvars -except CT_Corrected* gradients*

% Derive CT
CT_CorrectedArrayHigh = CT_CorrectedArray(CT_Corrected.Status == 'High', :);

% Compute and Fisher-transform CT covariance
CT_CovHigh = atanh(corr(CT_CorrectedArrayHigh));

% Nullify 'Inf' (along the diagonal)
CT_CovHigh(isinf(CT_CovHigh)) = 0;

% Create a gradient backbone using DME, normalised angle, and Procrustes 
% alignment
gradientBone = GradientMaps('approach', 'dm', 'kernel', 'na', 'alignment', 'Procrustes');

% Fit the matrix to the backbone with the total-sample gradients as the
% reference
gradientsCT_High = gradientBone.fit(CT_CovHigh, 'reference', gradientsArrayCT_Total);

% Derive the first 10 gradients
gradientsArrayCT_High = gradientsCT_High.aligned{1};

% Derive G1, aligned
G1_CT_High = gradientsCT_High.aligned{1}(:, 1);

save gradientsCT_High
load gradientsCT_Low

% Plot G1, aligned
G1_CT_HighConte69 = parcel_to_surface(G1_CT_High, 'glasser_360_conte69_fs_LR');
G1_CT_HighConte69(G1_CT_HighConte69 == 0) = (min(G1_CT_Low) + max(G1_CT_Low)) / 2;
figure('units', 'normalized', 'outerposition', [0.125 0 0.693 1])
title('CT G1, aligned: High (Low range)', 'FontSize', 30)
axis off
plot_cortical(G1_CT_HighConte69, 'surface_name', 'conte69', 'label_text', 'Loading', ...
     'color_range', [min(G1_CT_Low) max(G1_CT_Low)]);
colormap(flipud(brewermap([], 'PuOr')))
set(gca, 'FontSize', 25)

colourbars = findall(gcf, 'type', 'colorbar');
for cb = colourbars'
    set(cb, 'Ticks', [min(G1_CT_Low) max(G1_CT_Low)]);
    set(cb, 'TickLabels', {sprintf('%.3f', min(G1_CT_Low)) sprintf('%.3f', max(G1_CT_Low))});
end

%% 3.1 Correct SA for covariates: HCP

clc; clear; close all

load dataHCP
load data

clearvars -except dataHCP labelsGlasser*

dataHCP_Temp = dataHCP(:, {'ID' 'Age' 'IQ' 'eTIV'});
SA_CorrectedHCP = dataHCP(:, labelsGlasserSA);
SA_CorrectedHCP = [dataHCP_Temp SA_CorrectedHCP];
SA_CorrectedHCP_Array = nan(height(SA_CorrectedHCP), width(SA_CorrectedHCP) - width(dataHCP_Temp));
clearvars dataHCP_Temp

measure = 'Area';
covariates = '~ Age + IQ + eTIV';

measureOfInterest = strfind(SA_CorrectedHCP.Properties.VariableNames, measure);
measureOfInterestID = find(~cellfun(@isempty, measureOfInterest));

for parcelID = 1 : numel(measureOfInterestID)

    parcel = SA_CorrectedHCP.Properties.VariableNames{measureOfInterestID(parcelID)};

    fprintf('Correcting: %s\n', parcel)

    modelToBeFitted = sprintf('%s %s', parcel, covariates);
    model = fitlm(SA_CorrectedHCP, modelToBeFitted, 'RobustOpts', 'on');
    SA_CorrectedHCP_Array(:, parcelID) = model.Residuals.Raw;

end

SA_CorrectedHCP(:, labelsGlasserSA) = array2table(SA_CorrectedHCP_Array);

clearvars -except data* SA* labels*

save SA_CorrectedHCP

%% 3.2 Compute and plot SA G1: HCP

clc; clear; close all

load SA_CorrectedHCP

clearvars -except SA_Corrected*

% Compute and Fisher-transform SA covariance
SA_CovHCP = atanh(corr(SA_CorrectedHCP_Array));

% Nullify 'Inf' (along the diagonal)
SA_CovHCP(isinf(SA_CovHCP)) = 0;

% Create a gradient backbone using DME and normalised angle
gradientBone = GradientMaps('approach', 'dm', 'kernel', 'na');

% Fit the matrix to the backbone
gradientsSA_HCP = gradientBone.fit(SA_CovHCP);

% Derive the first 10 gradients
gradientsArraySA_HCP = gradientsSA_HCP.gradients{1};

% Derive G1
G1_SA_HCP = gradientsArraySA_HCP(:, 1);

save gradientsSA_HCP

% Plot G1
G1_SA_HCP_Conte69 = parcel_to_surface(G1_SA_HCP, 'glasser_360_conte69_fs_LR');
G1_SA_HCP_Conte69(G1_SA_HCP_Conte69 == 0) = (min(G1_SA_HCP) + max(G1_SA_HCP)) / 2;
figure('units', 'normalized', 'outerposition', [0.125 0 0.693 1])
title('SA G1: HCP', 'FontSize', 30)
axis off
plot_cortical(G1_SA_HCP_Conte69, 'surface_name', 'conte69', 'label_text', 'Loading');
colormap(flipud(brewermap([], 'BrBg')))
set(gca, 'FontSize', 25)

colourbars = findall(gcf, 'type', 'colorbar');
for cb = colourbars'
    set(cb, 'Ticks', [min(G1_SA_HCP) max(G1_SA_HCP)]);
    set(cb, 'TickLabels', {sprintf('%.3f', min(G1_SA_HCP)) sprintf('%.3f', max(G1_SA_HCP))});
end

%% 4.1 Correct SA for covariates: Total

clc; clear; close all

load data

clearvars -except data labelsGlasser*

dataTemp = data(:, {'URSI' 'PCL_Total' 'Status' 'Age' 'IQ' 'eTIV'});
SA_Corrected = data(:, labelsGlasserSA);
SA_Corrected = [dataTemp SA_Corrected];
SA_CorrectedArray = nan(height(SA_Corrected), width(SA_Corrected) - width(dataTemp));
clearvars dataTemp

measure = 'Area';
covariates = '~ Age + IQ + eTIV';

measureOfInterest = strfind(SA_Corrected.Properties.VariableNames, measure);
measureOfInterestURSI = find(~cellfun(@isempty, measureOfInterest));

for parcelID = 1 : numel(measureOfInterestURSI)

parcel = SA_Corrected.Properties.VariableNames{measureOfInterestURSI(parcelID)};

    fprintf('Correcting: %s\n', parcel)

    modelToBeFitted = sprintf('%s %s', parcel, covariates);
    model = fitlm(SA_Corrected, modelToBeFitted, 'RobustOpts', 'on');
    SA_CorrectedArray(:, parcelID) = model.Residuals.Raw;

end

SA_Corrected(:, labelsGlasserSA) = array2table(SA_CorrectedArray);

clearvars -except data* SA* labels*

save SA_Corrected

%% 4.2 Compute and plot SA G1: Total

clc; clear; close all

load SA_Corrected

clearvars -except SA_Corrected*

% Compute and Fisher-transform SA covariance
SA_CovTotal = atanh(corr(SA_CorrectedArray));

% Nullify 'Inf' (along the diagonal)
SA_CovTotal(isinf(SA_CovTotal)) = 0;

% Create a gradient backbone using DME and normalised angle
gradientBone = GradientMaps('approach', 'dm', 'kernel', 'na');

% Fit the matrix to the backbone
gradientsSA_Total = gradientBone.fit(SA_CovTotal);

% Derive the first 10 gradients
gradientsArraySA_Total = gradientsSA_Total.gradients{1};

% Derive G1
G1_SA_Total = gradientsArraySA_Total(:, 1);

save gradientsSA_Total

% Plot G1
G1_SA_TotalConte69 = parcel_to_surface(G1_SA_Total, 'glasser_360_conte69_fs_LR');
G1_SA_TotalConte69(G1_SA_TotalConte69 == 0) = (min(G1_SA_Total) + max(G1_SA_Total)) / 2;
figure('units', 'normalized', 'outerposition', [0.125 0 0.693 1])
title('SA G1: Total', 'FontSize', 30)
axis off
plot_cortical(G1_SA_TotalConte69, 'surface_name', 'conte69', 'label_text', 'Loading');
colormap(flipud(brewermap([], 'BrBg')))
set(gca, 'FontSize', 25)

colourbars = findall(gcf, 'type', 'colorbar');
for cb = colourbars'
    set(cb, 'Ticks', [min(G1_SA_Total) max(G1_SA_Total)]);
    set(cb, 'TickLabels', {sprintf('%.3f', min(G1_SA_Total)) sprintf('%.3f', max(G1_SA_Total))});
end

%% 4.3 Compute and plot SA G1: Low

clc; clear; close all

load SA_Corrected
load gradientsSA_Total

clearvars -except SA_Corrected* gradients*

% Derive SA
SA_CorrectedArrayLow = SA_CorrectedArray(SA_Corrected.Status == 'Low', :);

% Compute and Fisher-transform SA covariance
SA_CovLow = atanh(corr(SA_CorrectedArrayLow));

% Nullify 'Inf' (along the diagonal)
SA_CovLow(isinf(SA_CovLow)) = 0;

% Create a gradient backbone using DME, normalised angle, and Procrustes 
% alignment
gradientBone = GradientMaps('approach', 'dm', 'kernel', 'na', 'alignment', 'Procrustes');

% Fit the matrix to the backbone with the total-sample gradients as the
% reference
gradientsSA_Low = gradientBone.fit(SA_CovLow, 'reference', gradientsArraySA_Total);

% Derive the first 10 gradients
gradientsArraySA_Low = gradientsSA_Low.aligned{1};

% Derive G1, aligned
G1_SA_Low = gradientsSA_Low.aligned{1}(:, 1);

save gradientsSA_Low

% Plot G1, aligned
G1_SA_LowConte69 = parcel_to_surface(G1_SA_Low, 'glasser_360_conte69_fs_LR');
G1_SA_LowConte69(G1_SA_LowConte69 == 0) = (min(G1_SA_Low) + max(G1_SA_Low)) / 2;
figure('units', 'normalized', 'outerposition', [0.125 0 0.693 1])
title('SA G1, aligned: Low', 'FontSize', 30)
axis off
plot_cortical(G1_SA_LowConte69, 'surface_name', 'conte69', 'label_text', 'Loading');
colormap(flipud(brewermap([], 'BrBg')))
set(gca, 'FontSize', 25)

colourbars = findall(gcf, 'type', 'colorbar');
for cb = colourbars'
    set(cb, 'Ticks', [min(G1_SA_Low) max(G1_SA_Low)]);
    set(cb, 'TickLabels', {sprintf('%.3f', min(G1_SA_Low)) sprintf('%.3f', max(G1_SA_Low))});
end

%% 4.4 Compute and plot SA G1: High

clc; clear; close all

load SA_Corrected
load gradientsSA_Total

clearvars -except SA_Corrected* gradients*

% Derive SA
SA_CorrectedArrayHigh = SA_CorrectedArray(SA_Corrected.Status == 'High', :);

% Compute and Fisher-transform SA covariance
SA_CovHigh = atanh(corr(SA_CorrectedArrayHigh));

% Nullify 'Inf' (along the diagonal)
SA_CovHigh(isinf(SA_CovHigh)) = 0;

% Create a gradient backbone using DME, normalised angle, and Procrustes 
% alignment
gradientBone = GradientMaps('approach', 'dm', 'kernel', 'na', 'alignment', 'Procrustes');

% Fit the matrix to the backbone with the total-sample gradients as the
% reference
gradientsSA_High = gradientBone.fit(SA_CovHigh, 'reference', gradientsArraySA_Total);

% Derive the first 10 gradients
gradientsArraySA_High = gradientsSA_High.aligned{1};

% Derive G1, aligned
G1_SA_High = gradientsSA_High.aligned{1}(:, 1);

save gradientsSA_High
load gradientsSA_Low

% Plot G1, aligned
G1_SA_HighConte69 = parcel_to_surface(G1_SA_High, 'glasser_360_conte69_fs_LR');
G1_SA_HighConte69(G1_SA_HighConte69 == 0) = (min(G1_SA_Low) + max(G1_SA_Low)) / 2;
figure('units', 'normalized', 'outerposition', [0.125 0 0.693 1])
title('SA G1, aligned: High (Low range)', 'FontSize', 30)
axis off
plot_cortical(G1_SA_HighConte69, 'surface_name', 'conte69', 'label_text', 'Loading', ...
     'color_range', [min(G1_SA_Low) max(G1_SA_Low)]);
colormap(flipud(brewermap([], 'BrBg')))
set(gca, 'FontSize', 25)

colourbars = findall(gcf, 'type', 'colorbar');
for cb = colourbars'
    set(cb, 'Ticks', [min(G1_SA_Low) max(G1_SA_Low)]);
    set(cb, 'TickLabels', {sprintf('%.3f', min(G1_SA_Low)) sprintf('%.3f', max(G1_SA_Low))});
end

%% 5. HCP correlations and K-S tests

clc; clear; close all

load gradientsCT_HCP
load gradientsCT_Total
load gradientsSA_HCP
load gradientsSA_Total

rhoCT = corr(G1_CT_HCP, G1_CT_Total, 'Type', 'Spearman')
rhoSA = corr(G1_SA_HCP, G1_SA_Total, 'Type', 'Spearman')

rng default
spinP_CT = spin_test(G1_CT_HCP, G1_CT_Total, ...
     'parcellation_name', 'glasser_360', 'n_rot', 1000, 'Type', 'Spearman')

rng default
spinP_SA = spin_test(G1_SA_HCP, G1_SA_Total, ...
     'parcellation_name', 'glasser_360', 'n_rot', 1000, 'Type', 'Spearman')

save rhoCT_SA_HCP

clc; clear; close all

load gradientsCT_Low
load gradientsCT_High
load gradientsSA_Low
load gradientsSA_High

[~, KS_P_CT, KS_StatCT] = kstest2(G1_CT_High, G1_CT_Low)
[~, KS_P_SA, KS_StatSA] = kstest2(G1_SA_High, G1_SA_Low)