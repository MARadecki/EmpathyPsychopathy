%% Figure 4

%% 1.1 SA ~ Status

clc; clear; close all

load data

data = dataStatus;

clearvars -except data*

measure = '_Area';
IV = '~ Status';
covariates = ' + Age + IQ + eTIV';
covariatesN = count(covariates, '+') + 1;

measureOfInterest = strfind(data.Properties.VariableNames, measure);
measureOfInterestID = find(~cellfun(@isempty, measureOfInterest));

for parcelID = 1 : numel(measureOfInterestID)
    
    parcel = data.Properties.VariableNames{measureOfInterestID(parcelID)};
    data(:, {parcel}) = array2table(zscore(table2array(data(:, {parcel}))), 'VariableNames', {parcel});
    fprintf('Fitting: %s\n', parcel)

    modelToBeFitted = sprintf('%s %s', parcel, [IV covariates]);
    model = fitlm(data, modelToBeFitted, 'RobustOpts', 'on');
    modelParcelID = strcmp(model.CoefficientNames, 'Status_High');

    betasStatusSA(parcelID) = model.Coefficients.Estimate(modelParcelID);
    
    CI_Temporary = model.coefCI;
    confLowerStatusSA(parcelID) = CI_Temporary(modelParcelID, 1);
    confUpperStatusSA(parcelID) = CI_Temporary(modelParcelID, 2);
    
    P_ValuesStatusSA(parcelID) = model.Coefficients.pValue(modelParcelID);
    
    R_SquaredStatusSA(parcelID) = model.Rsquared.Adjusted;
    
end

P_ValuesStatusSA_FDR = mafdr(P_ValuesStatusSA, 'BHFDR', true);
P_ValuesStatusSA_FDR_Thr = P_ValuesStatusSA_FDR < 0.05;

save StatusSA

%% 1.2 SA ~ Status: Plot

clc; clear; close all

load StatusSA

range = [-0.3 0.3];

betasStatusSA_FDR = betasStatusSA;
betasStatusSA_FDR(~P_ValuesStatusSA_FDR_Thr) = 0;
betasStatusSA_FDR_Surface = parcel_to_surface(betasStatusSA_FDR, 'glasser_360_conte69_fs_LR');

figure('units', 'normalized', 'outerposition', [0.125 0 0.693 1])
T = tiledlayout(1, 4);
title(T, 'SA ~ Status', 'FontSize', 30, 'FontWeight', 'Bold', 'Color', 'Black')
axis off
plot_cortical(betasStatusSA_FDR_Surface, 'surface_name', 'conte69', 'label_text', 'β_Z', ...
'color_range', range, 'cmap', 'RdBu_r');
set(gca, 'FontSize', 20)

betasStatusSA_Surface = parcel_to_surface(betasStatusSA, 'glasser_360_conte69_fs_LR');

figure('units', 'normalized', 'outerposition', [0.125 0 0.693 1])
title({'SA: High vs low' sprintf('N = %.0f vs %.0f', sum(data.Status == 'High'), sum(data.Status == 'Low'))}, 'FontSize', 25)
axis off
plot_cortical(betasStatusSA_Surface, 'surface_name', 'conte69', 'label_text', 'β_Z', ...
'color_range', range, 'cmap', 'RdBu_r');
set(gca, 'FontSize', 20)

%% 2.1 Beta by class: Status

clc; clear; close all

load classesMesulamGlasser
load StatusSA

classNames = {'P' 'H' 'U' 'I'};

N_Groups = max(MesulamGlasser);

N_Tests = N_Groups * (N_Groups - 1) / 2;

IDX_Pair = 1;
for i = 1 : N_Groups
     for j = i + 1 : N_Groups
        
     map1 = betasStatusSA(MesulamGlasser == i);
     map2 = betasStatusSA(MesulamGlasser == j);

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

save betasStatusSA_Mesulam

%% 2.2 Beta by network: Status

clc; clear; close all

load networksYeoGlasser
load StatusSA

networkNames = {'V' 'S' 'DA' 'VA' 'L' 'F' 'DM'};

N_Groups = max(YeoGlasser);

N_Tests = N_Groups * (N_Groups - 1) / 2;

IDX_Pair = 1;
for i = 1 : N_Groups
     for j = i + 1 : N_Groups
        
     map1 = betasStatusSA(YeoGlasser == i);
     map2 = betasStatusSA(YeoGlasser == j);

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

save betasStatusSA_Yeo

%% 3.1 FDR-cluster overlap

clc; clear; close all

load StatusSA
load clustersSchurzGlasser

CognitiveOverlap = P_ValuesStatusSA_FDR_Thr' & logical(Cognitive);
CognitiveOverlapPerc = sum(CognitiveOverlap) / sum(P_ValuesStatusSA_FDR_Thr);

CogAffOverlap = P_ValuesStatusSA_FDR_Thr' & logical(CogAff);
CogAffOverlapPerc = sum(CogAffOverlap) / sum(P_ValuesStatusSA_FDR_Thr);

CogUniqueOverlap = P_ValuesStatusSA_FDR_Thr' & logical(CogUnique);
CogUniqueOverlapPerc = sum(CogUniqueOverlap) / sum(P_ValuesStatusSA_FDR_Thr);

AffUniqueOverlap = P_ValuesStatusSA_FDR_Thr' & logical(AffUnique);
AffUniqueOverlapPerc = sum(AffUniqueOverlap) / sum(P_ValuesStatusSA_FDR_Thr);

AffCogOverlap = P_ValuesStatusSA_FDR_Thr' & logical(AffCog);
AffCogOverlapPerc = sum(AffCogOverlap) / sum(P_ValuesStatusSA_FDR_Thr);

AffectiveOverlap = P_ValuesStatusSA_FDR_Thr' & logical(Affective);
AffectiveOverlapPerc = sum(AffectiveOverlap) / sum(P_ValuesStatusSA_FDR_Thr);

%

figure('units', 'normalized', 'outerposition', [0 0 1 1])
T = tiledlayout(9, 4);
title(T, 'Title', 'FontSize', 35, 'FontWeight', 'Bold', 'Color', 'Black')
subtitle(T, 'Subtitle', 'FontSize', 25, 'Color', 'Black')

nexttile([2 2])

hold on
line([-10 -10], [-10 -10], 'LineWidth', 1.5, 'Color', [0.25 0.25 0.25], 'LineStyle', '--')

barh(6, -CognitiveOverlapPerc, 'EdgeColor', 'b', 'LineWidth', 3, 'FaceColor', 'w', 'BarWidth', 0.45)
hold on
barh(5, -CogAffOverlapPerc, 'EdgeColor', 'b', 'LineWidth', 3, 'FaceColor', 'y', 'FaceAlpha', 1, 'BarWidth', 0.45)
hold on
barh(4, -CogUniqueOverlapPerc, 'EdgeColor', 'b', 'LineWidth', 3, 'FaceColor', 'g', 'FaceAlpha', 1, 'BarWidth', 0.45)
hold on
barh(3, -AffUniqueOverlapPerc, 'EdgeColor', 'r', 'LineWidth', 3, 'FaceColor', 'g', 'FaceAlpha', 1, 'BarWidth', 0.45)
hold on
barh(2, -AffCogOverlapPerc, 'EdgeColor', 'r', 'LineWidth', 3, 'FaceColor', 'y', 'FaceAlpha', 1, 'BarWidth', 0.45)
hold on
barh(1, -AffectiveOverlapPerc, 'EdgeColor', 'r', 'LineWidth', 3, 'FaceColor', 'w', 'BarWidth', 0.45)

set(gca, 'FontSize', 18)
xlabel('Scaled overlap with SA ~ Group at P_F_D_R < 0.05 [%]')
xlim([-AffectiveOverlapPerc 0])
xticks([-AffectiveOverlapPerc 0])
xticklabels({sprintf('%.0f', AffectiveOverlapPerc * 100) 0})
ylim([0 7])
yticks([1 : 6])
yticklabels(flip({'Cognitive' 'Cog: Pref' 'Cog: Unique' 'Aff: Unique' 'Aff: Pref' 'Affective'}))
box off

text(-CognitiveOverlapPerc - 0.0125, 6, sprintf('%.0f', CognitiveOverlapPerc*100), 'FontSize', 15, 'HorizontalAlignment', 'center', 'Color', 'b')
text(-CogAffOverlapPerc - 0.0125, 5, sprintf('%.0f', CogAffOverlapPerc*100), 'FontSize', 15, 'HorizontalAlignment', 'center', 'Color', 'b')
text(-CogUniqueOverlapPerc - 0.0125, 4, sprintf('%.0f', CogUniqueOverlapPerc*100), 'FontSize', 15, 'HorizontalAlignment', 'center', 'Color', 'b')
text(-AffUniqueOverlapPerc - 0.0125, 3, sprintf('%.0f', AffUniqueOverlapPerc*100), 'FontSize', 15, 'HorizontalAlignment', 'center', 'Color', 'r')
text(-AffCogOverlapPerc - 0.0125, 2, sprintf('%.0f', AffCogOverlapPerc*100), 'FontSize', 15, 'HorizontalAlignment', 'center', 'Color', 'r')

ax = gca;
ax.YAxisLocation = 'right';
ax.YAxis.FontSize = 15;

%% 3.2: Cluster-class/network overlap

clc; clear; close all

load StatusSA
load classesMesulamGlasser
load networksYeoGlasser
load clustersSchurzGlasser

CogAffOverlap = P_ValuesStatusSA_FDR_Thr' & logical(CogUnique);
CogAffOverlapPerc = sum(CogAffOverlap) / sum(P_ValuesStatusSA_FDR_Thr);

AffCogOverlap = P_ValuesStatusSA_FDR_Thr' & logical(AffUnique);
AffCogOverlapPerc = sum(AffCogOverlap) / sum(P_ValuesStatusSA_FDR_Thr);

%

MesulamLabels = {'P' 'H' 'U' 'I'};
YeoLabels = {'V' 'S' 'DA' 'VA' 'L' 'F' 'DM'};
MesulamRGB(1, :) = [];
YeoRGB(1, :) = [];

% Class
for i = 1 : max(MesulamGlasser)
     overlapCogAffClass(i) = sum(CogUnique & MesulamGlasser == i) / sum(MesulamGlasser == i);
     overlapAffCogClass(i) = sum(AffUnique & MesulamGlasser == i) / sum(MesulamGlasser == i);
end

% Network
for i = 1 : max(YeoGlasser)
     overlapCogAffNetwork(i) = sum(CogUnique & YeoGlasser == i) / sum(YeoGlasser == i);
     overlapAffCogNetwork(i) = sum(AffUnique & YeoGlasser == i) / sum(YeoGlasser == i);
end

%

figure('units', 'normalized', 'outerposition', [0 0 1 1])
T = tiledlayout(10, 4);

% Class

nexttile([4 1])

for i = 1 : 4
     hold on
     bar(i-0.125, overlapCogAffClass(i), 'EdgeColor', 'blue', 'LineWidth', 3, 'FaceColor', MesulamRGB(i, :), 'FaceAlpha', 0.75, 'BarWidth', 0.45)
     hold on
     bar(i+0.125, overlapAffCogClass(i), 'EdgeColor', 'red', 'LineWidth', 3, 'FaceColor', MesulamRGB(i, :), 'FaceAlpha', 0.75, 'BarWidth', 0.45)
end

bar(100, 100, 'EdgeColor', 'blue', 'LineWidth', 3, 'FaceColor', 'white')
bar(100, 100, 'EdgeColor', 'red', 'LineWidth', 3, 'FaceColor', 'white')

set(gca, 'FontSize', 18)
xlabel('Class')
xlim([0 5])
xticks([1 : 4])
xticklabels(MesulamLabels)
xtickangle(45)
ylim([0 max([overlapCogAffClass overlapAffCogClass])])
yticks([0 max([overlapCogAffClass overlapAffCogClass])])
yticklabels({0 sprintf('%.0f', max([overlapCogAffClass overlapAffCogClass]) * 100)})
box off

ax = gca;
ax.YAxisLocation = 'right';

% Network

nexttile([4 1])

for i = 1 : 7
     hold on
     bar(i-0.125, overlapCogAffNetwork(i), 'EdgeColor', 'blue', 'LineWidth', 3, 'FaceColor', YeoRGB(i, :), 'FaceAlpha', 0.75, 'BarWidth', 0.45)
     hold on
     bar(i+0.125, overlapAffCogNetwork(i), 'EdgeColor', 'red', 'LineWidth', 3, 'FaceColor', YeoRGB(i, :), 'FaceAlpha', 0.75, 'BarWidth', 0.45)
end

bar(100, 100, 'EdgeColor', 'blue', 'LineWidth', 3, 'FaceColor', 'white')
bar(100, 100, 'EdgeColor', 'red', 'LineWidth', 3, 'FaceColor', 'white')

set(gca, 'FontSize', 18)
xlabel('Network')
xlim([0 8])
xticks([1 : 7])
xticklabels(YeoLabels)
xtickangle(45)
ylabel({'Scaled overlap' 'with Cog: Unique & Aff: Unique [%]'})
ylim([0 max([overlapCogAffNetwork overlapAffCogNetwork])])
yticks([0 max([overlapCogAffNetwork overlapAffCogNetwork])])
yticklabels({0 sprintf('%.0f', max([overlapCogAffNetwork overlapAffCogNetwork]) * 100)})
box off

ax = gca;
ax.YAxisLocation = 'right';

%% 3.3 FDR-Neurosynth overlap

clc; clear; close all

load StatusSA
load NS_24

for i = 1 : height(NS_Labels)

    NS_Label = NS_Labels{i};
    NS_Data = eval(NS_Label);
    NS_Overlap(i) = sum(P_ValuesStatusSA_FDR_Thr' == 1 & logical(NS_Data) == 1) / sum(P_ValuesStatusSA_FDR_Thr);

end

tableStatusNS = array2table(NS_Overlap', 'VariableNames', {'NS_Overlap'});
tableStatusNS.NS_Labels = NS_Labels;
tableStatusNS.NS_Labels = replace(tableStatusNS.NS_Labels, '_', ' ');
tableStatusNS = sortrows(tableStatusNS, 'NS_Overlap', 'ascend');
tableStatusNS.NS_Labels = replace(tableStatusNS.NS_Labels, 'autobiographical memory', 'autobio. memory');

%

figure('units', 'normalized', 'outerposition', [0 0 1 1])
T = tiledlayout(10, 1);

nexttile([4 1])

line([-50 50], [100 100], 'LineWidth', 1.5, 'Color', [0.85 0.85 0.85], 'LineStyle', '--')
hold on

barNetwork = bar(-tableStatusNS.NS_Overlap, 'FaceColor', [0.85 0.85 0.85], 'EdgeColor', 'black', 'LineWidth', 3, 'BarWidth', 0.65);
set(gca, 'FontSize', 18)
xlabel('Neurosynth term')
xlim([0 height(tableStatusNS)+1])
xticks(1 : height(tableStatusNS))
xticklabels(tableStatusNS.NS_Labels)
xtickangle(35)
ylabel({'Scaled overlap' 'with P_F_D_R < 0.05 [%]'})
ylim([-max(tableStatusNS.NS_Overlap) 0])
yticks([-max(tableStatusNS.NS_Overlap) 0])
yticklabels({sprintf('%.0f', max(tableStatusNS.NS_Overlap) * 100) 0})
box off

ax = gca;
ax.TickLength = [ax.TickLength(1)/1.75 ax.TickLength(1)/1.75];
ax.YAxisLocation = 'right';

colourMap = bone(height(tableStatusNS));
set(barNetwork, 'FaceColor', 'flat');
set(barNetwork, 'CData', colourMap);