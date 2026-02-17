%% Figure 2

%% 1.1 IRI-PT ~ SA

clc; clear; close all

load data

clearvars -except data*

data.IRI_PT = zscore(data.IRI_PT);

DV = 'IRI_PT';
measure = '_Area';
covariates = '~ Age + IQ + eTIV + ';

measureOfInterest = strfind(data.Properties.VariableNames, measure);
measureOfInterestID = find(~cellfun(@isempty, measureOfInterest));

for parcelID = 1 : numel(measureOfInterestID)
    
    parcel = data.Properties.VariableNames{measureOfInterestID(parcelID)};
    data(:, {parcel}) = array2table(zscore(table2array(data(:, {parcel}))), 'VariableNames', {parcel});
    fprintf('Fitting: %s\n', parcel)

    modelToBeFitted = sprintf('%s %s', DV, [covariates parcel]);
    model = fitlm(data, modelToBeFitted, 'RobustOpts', 'on');
    modelParcelID = strcmp(model.CoefficientNames, parcel);

    betasIRI_PT_SA(parcelID) = model.Coefficients.Estimate(modelParcelID);
    
    CI_Temporary = model.coefCI;
    confLowerIRI_PT_SA(parcelID) = CI_Temporary(modelParcelID, 1);
    confUpperIRI_PT_SA(parcelID) = CI_Temporary(modelParcelID, 2);
    
    P_ValuesIRI_PT_SA(parcelID) = model.Coefficients.pValue(modelParcelID);
    
    R_SquaredIRI_PT_SA(parcelID) = model.Rsquared.Adjusted;
    
end

P_ValuesIRI_PT_SA_FDR = mafdr(P_ValuesIRI_PT_SA, 'BHFDR', true);
P_ValuesIRI_PT_SA_FDR_Thr = P_ValuesIRI_PT_SA_FDR < 0.05;

clearvars -except data DV covariates betas* conf* P_Values* R_Squared* *Temporary

save IRI_PT_SA

%% 1.2 IRI-EC ~ SA

clc; clear; close all

load data

clearvars -except data*

data.IRI_EC = zscore(data.IRI_EC);

DV = 'IRI_EC';
measure = '_Area';
covariates = '~ Age + IQ + eTIV + ';

measureOfInterest = strfind(data.Properties.VariableNames, measure);
measureOfInterestID = find(~cellfun(@isempty, measureOfInterest));

for parcelID = 1 : numel(measureOfInterestID)
    
    parcel = data.Properties.VariableNames{measureOfInterestID(parcelID)};
    data(:, {parcel}) = array2table(zscore(table2array(data(:, {parcel}))), 'VariableNames', {parcel});
    fprintf('Fitting: %s\n', parcel)

    modelToBeFitted = sprintf('%s %s', DV, [covariates parcel]);
    model = fitlm(data, modelToBeFitted, 'RobustOpts', 'on');
    modelParcelID = strcmp(model.CoefficientNames, parcel);

    betasIRI_EC_SA(parcelID) = model.Coefficients.Estimate(modelParcelID);
    
    CI_Temporary = model.coefCI;
    confLowerIRI_EC_SA(parcelID) = CI_Temporary(modelParcelID, 1);
    confUpperIRI_EC_SA(parcelID) = CI_Temporary(modelParcelID, 2);
    
    P_ValuesIRI_EC_SA(parcelID) = model.Coefficients.pValue(modelParcelID);
    
    R_SquaredIRI_EC_SA(parcelID) = model.Rsquared.Adjusted;
    
end

P_ValuesIRI_EC_SA_FDR = mafdr(P_ValuesIRI_EC_SA, 'BHFDR', true);
P_ValuesIRI_EC_SA_FDR_Thr = P_ValuesIRI_EC_SA_FDR < 0.05;

clearvars -except data DV covariates betas* conf* P_Values* R_Squared* *Temporary

save IRI_EC_SA

%% 1.3 PCL-R F1 ~ SA

clc; clear; close all

load data

clearvars -except data*

data.PCL_Factor1 = zscore(data.PCL_Factor1);

DV = 'PCL_Factor1';
measure = '_Area';
covariates = '~ Age + IQ + eTIV + ';

measureOfInterest = strfind(data.Properties.VariableNames, measure);
measureOfInterestID = find(~cellfun(@isempty, measureOfInterest));

for parcelID = 1 : numel(measureOfInterestID)
    
    parcel = data.Properties.VariableNames{measureOfInterestID(parcelID)};
    data(:, {parcel}) = array2table(zscore(table2array(data(:, {parcel}))), 'VariableNames', {parcel});
    fprintf('Fitting: %s\n', parcel)

    modelToBeFitted = sprintf('%s %s', DV, [covariates parcel]);
    model = fitlm(data, modelToBeFitted, 'RobustOpts', 'on');
    modelParcelID = strcmp(model.CoefficientNames, parcel);

    betasPCL_Factor1_SA(parcelID) = model.Coefficients.Estimate(modelParcelID);
    
    CI_Temporary = model.coefCI;
    confLowerPCL_Factor1_SA(parcelID) = CI_Temporary(modelParcelID, 1);
    confUpperPCL_Factor1_SA(parcelID) = CI_Temporary(modelParcelID, 2);
    
    P_ValuesPCL_Factor1_SA(parcelID) = model.Coefficients.pValue(modelParcelID);
    
    R_SquaredPCL_Factor1_SA(parcelID) = model.Rsquared.Adjusted;
    
end

P_ValuesPCL_Factor1_SA_FDR = mafdr(P_ValuesPCL_Factor1_SA, 'BHFDR', true);
P_ValuesPCL_Factor1_SA_FDR_Thr = P_ValuesPCL_Factor1_SA_FDR < 0.05;

clearvars -except data DV covariates betas* conf* P_Values* R_Squared* *Temporary

save PCL_Factor1_SA

%% 1.4 PCL-R F2 ~ SA

clc; clear; close all

load data

clearvars -except data*

data = rmmissing(data, 'DataVariables', 'PCL_Factor2');
data.PCL_Factor2 = zscore(data.PCL_Factor2);

DV = 'PCL_Factor2';
measure = '_Area';
covariates = '~ Age + IQ + eTIV + ';

measureOfInterest = strfind(data.Properties.VariableNames, measure);
measureOfInterestID = find(~cellfun(@isempty, measureOfInterest));

for parcelID = 1 : numel(measureOfInterestID)
    
    parcel = data.Properties.VariableNames{measureOfInterestID(parcelID)};
    data(:, {parcel}) = array2table(zscore(table2array(data(:, {parcel}))), 'VariableNames', {parcel});
    fprintf('Fitting: %s\n', parcel)

    modelToBeFitted = sprintf('%s %s', DV, [covariates parcel]);
    model = fitlm(data, modelToBeFitted, 'RobustOpts', 'on');
    modelParcelID = strcmp(model.CoefficientNames, parcel);

    betasPCL_Factor2_SA(parcelID) = model.Coefficients.Estimate(modelParcelID);
    
    CI_Temporary = model.coefCI;
    confLowerPCL_Factor2_SA(parcelID) = CI_Temporary(modelParcelID, 1);
    confUpperPCL_Factor2_SA(parcelID) = CI_Temporary(modelParcelID, 2);
    
    P_ValuesPCL_Factor2_SA(parcelID) = model.Coefficients.pValue(modelParcelID);
    
    R_SquaredPCL_Factor2_SA(parcelID) = model.Rsquared.Adjusted;
    
end

P_ValuesPCL_Factor2_SA_FDR = mafdr(P_ValuesPCL_Factor2_SA, 'BHFDR', true);
P_ValuesPCL_Factor2_SA_FDR_Thr = P_ValuesPCL_Factor2_SA_FDR < 0.05;

clearvars -except data DV covariates betas* conf* P_Values* R_Squared* *Temporary

save PCL_Factor2_SA

%% 2.1 IRI-PT ~ SA: Plot

clc; clear; close all

load IRI_PT_SA

betasIRI_PT_SA_FDR = betasIRI_PT_SA;
betasIRI_PT_SA_FDR(~P_ValuesIRI_PT_SA_FDR_Thr) = 0;
betasIRI_PT_SA_FDR_Surface = parcel_to_surface(betasIRI_PT_SA_FDR, 'glasser_360_conte69_fs_LR');

figure('units', 'normalized', 'outerposition', [0.125 0 0.693 1])
T = tiledlayout(1, 4);
title(T, sprintf('%s ~ SA (FDR)', DV), 'FontSize', 30, 'FontWeight', 'Bold', 'Color', 'Black')
axis off
plot_cortical(betasIRI_PT_SA_FDR_Surface, 'surface_name', 'conte69', 'label_text', 'β_Z', ...
'color_range', [-0.15 0.15], 'cmap', 'RdBu_r');
set(gca, 'FontSize', 20)

IRI_PT_SA_Surface = parcel_to_surface(betasIRI_PT_SA, 'glasser_360_conte69_fs_LR');

figure('units', 'normalized', 'outerposition', [0.125 0 0.693 1])
T = tiledlayout(1, 4);
title(T, sprintf('%s ~ SA', DV), 'FontSize', 30, 'FontWeight', 'Bold', 'Color', 'Black')
axis off
plot_cortical(IRI_PT_SA_Surface, 'surface_name', 'conte69', 'label_text', 'β_Z', ...
'color_range', [-0.15 0.15], 'cmap', 'RdBu_r');
set(gca, 'FontSize', 20)

%% 2.2 IRI-EC ~ SA: Plot

clc; clear; close all

load IRI_EC_SA

betasIRI_EC_SA_FDR = betasIRI_EC_SA;
betasIRI_EC_SA_FDR(~P_ValuesIRI_EC_SA_FDR_Thr) = 0;
betasIRI_EC_SA_FDR_Surface = parcel_to_surface(betasIRI_EC_SA_FDR, 'glasser_360_conte69_fs_LR');

figure('units', 'normalized', 'outerposition', [0.125 0 0.693 1])
T = tiledlayout(1, 4);
title(T, sprintf('%s ~ SA (FDR)', DV), 'FontSize', 30, 'FontWeight', 'Bold', 'Color', 'Black')
axis off
plot_cortical(betasIRI_EC_SA_FDR_Surface, 'surface_name', 'conte69', 'label_text', 'β_Z', ...
'color_range', [-0.15 0.15], 'cmap', 'RdBu_r');
set(gca, 'FontSize', 20)

IRI_EC_SA_Surface = parcel_to_surface(betasIRI_EC_SA, 'glasser_360_conte69_fs_LR');

figure('units', 'normalized', 'outerposition', [0.125 0 0.693 1])
T = tiledlayout(1, 4);
title(T, sprintf('%s ~ SA', DV), 'FontSize', 30, 'FontWeight', 'Bold', 'Color', 'Black')
axis off
plot_cortical(IRI_EC_SA_Surface, 'surface_name', 'conte69', 'label_text', 'β_Z', ...
'color_range', [-0.15 0.15], 'cmap', 'RdBu_r');
set(gca, 'FontSize', 20)

%% 2.3 PCL-R F1 ~ SA: Plot

clc; clear; close all

load PCL_Factor1_SA

betasPCL_Factor1_SA_FDR = betasPCL_Factor1_SA;
betasPCL_Factor1_SA_FDR(~P_ValuesPCL_Factor1_SA_FDR_Thr) = 0;
betasPCL_Factor1_SA_FDR_Surface = parcel_to_surface(betasPCL_Factor1_SA_FDR, 'glasser_360_conte69_fs_LR');

figure('units', 'normalized', 'outerposition', [0.125 0 0.693 1])
T = tiledlayout(1, 4);
title(T, sprintf('%s ~ SA (FDR)', DV), 'FontSize', 30, 'FontWeight', 'Bold', 'Color', 'Black')
axis off
plot_cortical(betasPCL_Factor1_SA_FDR_Surface, 'surface_name', 'conte69', 'label_text', 'β_Z', ...
'color_range', [-0.15 0.15], 'cmap', 'RdBu_r');
set(gca, 'FontSize', 20)

PCL_Factor1_SA_Surface = parcel_to_surface(betasPCL_Factor1_SA, 'glasser_360_conte69_fs_LR');

figure('units', 'normalized', 'outerposition', [0.125 0 0.693 1])
T = tiledlayout(1, 4);
title(T, sprintf('%s ~ SA', DV), 'FontSize', 30, 'FontWeight', 'Bold', 'Color', 'Black')
axis off
plot_cortical(PCL_Factor1_SA_Surface, 'surface_name', 'conte69', 'label_text', 'β_Z', ...
'color_range', [-0.15 0.15], 'cmap', 'RdBu_r');
set(gca, 'FontSize', 20)

%% 2.4 PCL-R F2 ~ SA: Plot

clc; clear; close all

load PCL_Factor2_SA

betasPCL_Factor2_SA_FDR = betasPCL_Factor2_SA;
betasPCL_Factor2_SA_FDR(~P_ValuesPCL_Factor2_SA_FDR_Thr) = 0;
betasPCL_Factor2_SA_FDR_Surface = parcel_to_surface(betasPCL_Factor2_SA_FDR, 'glasser_360_conte69_fs_LR');

figure('units', 'normalized', 'outerposition', [0.125 0 0.693 1])
T = tiledlayout(1, 4);
title(T, sprintf('%s ~ SA (FDR)', DV), 'FontSize', 30, 'FontWeight', 'Bold', 'Color', 'Black')
axis off
plot_cortical(betasPCL_Factor2_SA_FDR_Surface, 'surface_name', 'conte69', 'label_text', 'β_Z', ...
'color_range', [-0.15 0.15], 'cmap', 'RdBu_r');
set(gca, 'FontSize', 20)

PCL_Factor2_SA_Surface = parcel_to_surface(betasPCL_Factor2_SA, 'glasser_360_conte69_fs_LR');

figure('units', 'normalized', 'outerposition', [0.125 0 0.693 1])
T = tiledlayout(1, 4);
title(T, sprintf('%s ~ SA', DV), 'FontSize', 30, 'FontWeight', 'Bold', 'Color', 'Black')
axis off
plot_cortical(PCL_Factor2_SA_Surface, 'surface_name', 'conte69', 'label_text', 'β_Z', ...
'color_range', [-0.15 0.15], 'cmap', 'RdBu_r');
set(gca, 'FontSize', 20)

%% 3.1 Beta by class: IRI-PT

clc; clear; close all

load classesMesulamGlasser
load IRI_PT_SA

classNames = {'P' 'H' 'U' 'I'};

N_Groups = max(MesulamGlasser);
N_Tests = N_Groups * (N_Groups - 1) / 2;

IDX_Pair = 1;
for i = 1 : N_Groups
     for j = i + 1 : N_Groups
     
     map1 = betasIRI_PT_SA(MesulamGlasser == i);
     map2 = betasIRI_PT_SA(MesulamGlasser == j);

     P_Values(IDX_Pair) = ranksum(map1, map2);

     comparisonNames{IDX_Pair, 1} = classNames{i};
     comparisonNames{IDX_Pair, 2} = classNames{j};

     IDX_Pair = IDX_Pair + 1;
     
     end
end

P_ValuesBon = min(1, P_Values * N_Tests);

disp("Significant comparisons following Bonferroni's correction:");
for k = 1 : N_Tests
    if P_ValuesBon(k) < 0.05
        fprintf('%s vs %s: P = %.4f\n', comparisonNames{k, 1}, comparisonNames{k, 2}, P_ValuesBon(k));
    end
end

save IRI_PT_SA_Mesulam

%% 3.2 Beta by class: IRI-EC

clc; clear; close all

load classesMesulamGlasser
load IRI_EC_SA

classNames = {'P' 'H' 'U' 'I'};

N_Groups = max(MesulamGlasser);
N_Tests = N_Groups * (N_Groups - 1) / 2;

IDX_Pair = 1;
for i = 1 : N_Groups
     for j = i + 1 : N_Groups
     
     map1 = betasIRI_EC_SA(MesulamGlasser == i);
     map2 = betasIRI_EC_SA(MesulamGlasser == j);

     P_Values(IDX_Pair) = ranksum(map1, map2);

     comparisonNames{IDX_Pair, 1} = classNames{i};
     comparisonNames{IDX_Pair, 2} = classNames{j};

     IDX_Pair = IDX_Pair + 1;
     
     end
end

P_ValuesBon = min(1, P_Values * N_Tests);

disp("Significant comparisons following Bonferroni's correction:");
for k = 1 : N_Tests
    if P_ValuesBon(k) < 0.05
        fprintf('%s vs %s: P = %.4f\n', comparisonNames{k, 1}, comparisonNames{k, 2}, P_ValuesBon(k));
    end
end

save IRI_EC_SA_Mesulam

%% 3.3 Beta by class: PCL-R F1

clc; clear; close all

load classesMesulamGlasser
load PCL_Factor1_SA

classNames = {'P' 'H' 'U' 'I'};

N_Groups = max(MesulamGlasser);
N_Tests = N_Groups * (N_Groups - 1) / 2;

IDX_Pair = 1;
for i = 1 : N_Groups
     for j = i + 1 : N_Groups
     
     map1 = betasPCL_Factor1_SA(MesulamGlasser == i);
     map2 = betasPCL_Factor1_SA(MesulamGlasser == j);

     P_Values(IDX_Pair) = ranksum(map1, map2);

     comparisonNames{IDX_Pair, 1} = classNames{i};
     comparisonNames{IDX_Pair, 2} = classNames{j};

     IDX_Pair = IDX_Pair + 1;
     
     end
end

P_ValuesBon = min(1, P_Values * N_Tests);

disp("Significant comparisons following Bonferroni's correction:");
for k = 1 : N_Tests
    if P_ValuesBon(k) < 0.05
        fprintf('%s vs %s: P = %.4f\n', comparisonNames{k, 1}, comparisonNames{k, 2}, P_ValuesBon(k));
    end
end

save PCL_Factor1_SA_Mesulam

%% 3.4 Beta by class: PCL-R F2

clc; clear; close all

load classesMesulamGlasser
load PCL_Factor2_SA

classNames = {'P' 'H' 'U' 'I'};

N_Groups = max(MesulamGlasser);
N_Tests = N_Groups * (N_Groups - 1) / 2;

IDX_Pair = 1;
for i = 1 : N_Groups
     for j = i + 1 : N_Groups
     
     map1 = betasPCL_Factor2_SA(MesulamGlasser == i);
     map2 = betasPCL_Factor2_SA(MesulamGlasser == j);

     P_Values(IDX_Pair) = ranksum(map1, map2);

     comparisonNames{IDX_Pair, 1} = classNames{i};
     comparisonNames{IDX_Pair, 2} = classNames{j};

     IDX_Pair = IDX_Pair + 1;
     
     end
end

P_ValuesBon = min(1, P_Values * N_Tests);

disp("Significant comparisons following Bonferroni's correction:");
for k = 1 : N_Tests
    if P_ValuesBon(k) < 0.05
        fprintf('%s vs %s: P = %.4f\n', comparisonNames{k, 1}, comparisonNames{k, 2}, P_ValuesBon(k));
    end
end

save PCL_Factor2_SA_Mesulam

%% 4.1 Beta by network: IRI-PT

clc; clear; close all

load networksYeoGlasser
load IRI_PT_SA

networkNames = {'V' 'S' 'DA' 'VA' 'L' 'F' 'DM'};

N_Groups = max(YeoGlasser);
N_Tests = N_Groups * (N_Groups - 1) / 2;

IDX_Pair = 1;
for i = 1 : N_Groups
     for j = i + 1 : N_Groups
     
     map1 = betasIRI_PT_SA(YeoGlasser == i);
     map2 = betasIRI_PT_SA(YeoGlasser == j);

     P_Values(IDX_Pair) = ranksum(map1, map2);

     comparisonNames{IDX_Pair, 1} = networkNames{i};
     comparisonNames{IDX_Pair, 2} = networkNames{j};

     IDX_Pair = IDX_Pair + 1;
     
     end
end

P_ValuesBon = min(1, P_Values * N_Tests);

disp("Significant comparisons following Bonferroni's correction:");
for k = 1 : N_Tests
    if P_ValuesBon(k) < 0.05
        fprintf('%s vs %s: P = %.4f\n', comparisonNames{k, 1}, comparisonNames{k, 2}, P_ValuesBon(k));
    end
end

save IRI_PT_SA_Yeo

%% 4.2 Beta by network: IRI-EC

clc; clear; close all

load networksYeoGlasser
load IRI_EC_SA

networkNames = {'V' 'S' 'DA' 'VA' 'L' 'F' 'DM'};

N_Groups = max(YeoGlasser);
N_Tests = N_Groups * (N_Groups - 1) / 2;

IDX_Pair = 1;
for i = 1 : N_Groups
     for j = i + 1 : N_Groups
     
     map1 = betasIRI_EC_SA(YeoGlasser == i);
     map2 = betasIRI_EC_SA(YeoGlasser == j);

     P_Values(IDX_Pair) = ranksum(map1, map2);

     comparisonNames{IDX_Pair, 1} = networkNames{i};
     comparisonNames{IDX_Pair, 2} = networkNames{j};

     IDX_Pair = IDX_Pair + 1;
        
     end
end

P_ValuesBon = min(1, P_Values * N_Tests);

disp("Significant comparisons following Bonferroni's correction:");
for k = 1 : N_Tests
    if P_ValuesBon(k) < 0.05
        fprintf('%s vs %s: P = %.4f\n', comparisonNames{k, 1}, comparisonNames{k, 2}, P_ValuesBon(k));
    end
end

save IRI_EC_SA_Yeo

%% 4.3 Beta by network: PCL-R F1

clc; clear; close all

load networksYeoGlasser
load PCL_Factor1_SA

networkNames = {'V' 'S' 'DA' 'VA' 'L' 'F' 'DM'};

N_Groups = max(YeoGlasser);
N_Tests = N_Groups * (N_Groups - 1) / 2;

IDX_Pair = 1;
for i = 1 : N_Groups
     for j = i + 1 : N_Groups
     
     map1 = betasPCL_Factor1_SA(YeoGlasser == i);
     map2 = betasPCL_Factor1_SA(YeoGlasser == j);

     P_Values(IDX_Pair) = ranksum(map1, map2);

     comparisonNames{IDX_Pair, 1} = networkNames{i};
     comparisonNames{IDX_Pair, 2} = networkNames{j};

     IDX_Pair = IDX_Pair + 1;
        
     end
end

P_ValuesBon = min(1, P_Values * N_Tests);

disp("Significant comparisons following Bonferroni's correction:");
for k = 1 : N_Tests
    if P_ValuesBon(k) < 0.05
        fprintf('%s vs %s: P = %.4f\n', comparisonNames{k, 1}, comparisonNames{k, 2}, P_ValuesBon(k));
    end
end

save PCL_Factor1_SA_Yeo

%% 4.4 Beta by network: PCL-R F2

clc; clear; close all

load networksYeoGlasser
load PCL_Factor2_SA

networkNames = {'V' 'S' 'DA' 'VA' 'L' 'F' 'DM'};

N_Groups = max(YeoGlasser);
N_Tests = N_Groups * (N_Groups - 1) / 2;

IDX_Pair = 1;
for i = 1 : N_Groups
     for j = i + 1 : N_Groups
     
     map1 = betasPCL_Factor2_SA(YeoGlasser == i);
     map2 = betasPCL_Factor2_SA(YeoGlasser == j);

     P_Values(IDX_Pair) = ranksum(map1, map2);

     comparisonNames{IDX_Pair, 1} = networkNames{i};
     comparisonNames{IDX_Pair, 2} = networkNames{j};

     IDX_Pair = IDX_Pair + 1;
        
     end
end

P_ValuesBon = min(1, P_Values * N_Tests);

disp("Significant comparisons following Bonferroni's correction:");
for k = 1 : N_Tests
    if P_ValuesBon(k) < 0.05
        fprintf('%s vs %s: P = %.4f\n', comparisonNames{k, 1}, comparisonNames{k, 2}, P_ValuesBon(k));
    end
end

save PCL_Factor2_SA_Yeo