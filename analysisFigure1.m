%% Figure 1

%% 1.1 IRI subscales & PCL-R factors

clc; clear; close all

load data

figure('units', 'normalized', 'outerposition', [0 0 1 1])
T = tiledlayout(2, 4);

% 1. IRI-PT

nexttile
H_IRI_PT = histogram(data.IRI_PT, 'BinWidth', 1, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'white');
H_IRI_PT.LineWidth = 1.5;
hold on
set(gca, 'FontSize', 18)
title('IRI-PT', 'FontSize', 25)
subtitle(sprintf('N = %.0f', sum(~isnan(data.IRI_PT))), 'FontSize', 20)
xlabel('Perspective Taking')
H_IRI_PT.BinLimits = [-0.5 28.5];
xlim([-2 30])
xticks([0 median(data.IRI_PT) 28])
ylabel('Frequency')
ylim([0 max(H_IRI_PT.BinCounts)])
yticks([0 max(H_IRI_PT.BinCounts)])
box off
axis square

line([median(data.IRI_PT) median(data.IRI_PT)], [0 100], 'LineWidth', 2, 'Color', 'black', 'LineStyle', '--')

[F_IRI_PT, XI_IRI_PT] = ksdensity(data.IRI_PT, 'BandWidth', 1.5);
plot(XI_IRI_PT, F_IRI_PT * sum(H_IRI_PT.BinWidth * H_IRI_PT.BinCounts), 'Color', [0.5 0.5 0.5], 'LineWidth', 4);

legendInfo = get(gca, 'Children');
legend(legendInfo(end - 1), 'Med', 'FontSize', 15, 'Location', 'NorthWest', 'box', 'off')

% 2. IRI-EC

nexttile
H_IRI_EC = histogram(data.IRI_EC, 'BinWidth', 1, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'white');
H_IRI_EC.LineWidth = 1.5;
hold on
set(gca, 'FontSize', 18)
title('IRI-EC', 'FontSize', 25)
subtitle(sprintf('N = %.0f', sum(~isnan(data.IRI_EC))), 'FontSize', 20)
xlabel('Empathic Concern')
H_IRI_EC.BinLimits = [-0.5 28.5];
xlim([-2 30])
xticks([0 median(data.IRI_EC) 28])
ylim([0 max(H_IRI_EC.BinCounts)])
yticks([0 max(H_IRI_EC.BinCounts)])
box off
axis square

line([median(data.IRI_EC) median(data.IRI_EC)], [0 100], 'LineWidth', 2, 'Color', 'black', 'LineStyle', '--')

[F_IRI_EC, XI_IRI_EC] = ksdensity(data.IRI_EC, 'BandWidth', 1.5);
plot(XI_IRI_EC, F_IRI_EC * sum(H_IRI_EC.BinWidth * H_IRI_EC.BinCounts), 'Color', [0.5 0.5 0.5], 'LineWidth', 4);

% 3. PCL-R F1

nexttile
H_PCL_Factor1 = histogram(data.PCL_Factor1, 'BinWidth', 1, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'white');
H_PCL_Factor1.LineWidth = 1.5;
hold on
set(gca, 'FontSize', 18)
title('PCL-R F1', 'FontSize', 25)
subtitle(sprintf('N = %.0f', sum(~isnan(data.PCL_Factor1))))
xlabel('Interpersonal/Affective')
H_PCL_Factor1.BinLimits = [-0.5 16.5];
xlim([-2 18])
xticks([0 median(data.PCL_Factor1) 16])
ylim([0 max(H_PCL_Factor1.BinCounts)])
yticks([0 max(H_PCL_Factor1.BinCounts)])
box off
axis square

line([median(data.PCL_Factor1) median(data.PCL_Factor1)], [0 100], 'LineWidth', 2, 'Color', 'black', 'LineStyle', '--')

[F_PCL_Factor1, XI_PCL_Factor1] = ksdensity(data.PCL_Factor1, 'BandWidth', 1.5);
plot(XI_PCL_Factor1, F_PCL_Factor1 * sum(H_PCL_Factor1.BinWidth * H_PCL_Factor1.BinCounts), 'Color', [0.5 0.5 0.5], 'LineWidth', 4);

% 4. PCL-R F2

nexttile
H_PCL_Factor2 = histogram(data.PCL_Factor2, 'BinWidth', 1, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'white');
H_PCL_Factor2.LineWidth = 1.5;
hold on
set(gca, 'FontSize', 18)
title('PCL-R F2', 'FontSize', 25)
subtitle(sprintf('N = %.0f', sum(~isnan(data.PCL_Factor2))), 'FontSize', 20)
xlabel('Lifestyle/Antisocial')
H_PCL_Factor2.BinLimits = [-0.5 20.5];
xlim([-2 22])
xticks([0 nanmedian(data.PCL_Factor2) 20])
ylim([0 max(H_PCL_Factor2.BinCounts)])
yticks([0 max(H_PCL_Factor2.BinCounts)])
box off
axis square

line([nanmedian(data.PCL_Factor2) nanmedian(data.PCL_Factor2)], [0 100], 'LineWidth', 2, 'Color', 'black', 'LineStyle', '--')

[F_PCL_Factor2, XI_PCL_Factor2] = ksdensity(data.PCL_Factor2, 'BandWidth', 1.5);
plot(XI_PCL_Factor2, F_PCL_Factor2 * sum(H_PCL_Factor2.BinWidth * H_PCL_Factor2.BinCounts), 'Color', [0.5 0.5 0.5], 'LineWidth', 4);

%% 1.2: IRI subscales ~ PCL-R factors

clc; clear; close all

load data

data.IRI_PT = zscore(data.IRI_PT);
data.IRI_EC = zscore(data.IRI_EC);
data.PCL_Factor1 = zscore(data.PCL_Factor1);

dataF2 = rmmissing(data, 'DataVariables', 'PCL_Factor2');

dataF2.IRI_PT = zscore(dataF2.IRI_PT);
dataF2.IRI_EC = zscore(dataF2.IRI_EC);
dataF2.PCL_Factor2 = zscore(dataF2.PCL_Factor2);

figure('units', 'normalized', 'outerposition', [0 0 1 1])
T = tiledlayout(2, 4);

% 1. IRI-PT ~ PCL-R F1

nexttile

line([0 0], [-3.5 3.5], 'LineWidth', 1.5, 'Color', [0.85 0.85 0.85], 'LineStyle', '--')
line([-3.5 3.5], [0 0], 'LineWidth', 1.5, 'Color', [0.85 0.85 0.85], 'LineStyle', '--')
hold on

rng default
spreadRng = 0.055 * randn(height(data), 1);

modelPT_F1 = fitlm(data, 'IRI_PT ~ PCL_Factor1 + Age + IQ', 'RobustOpts', 'on');
modelPT_F1_Plot = plotAdjustedResponse(modelPT_F1, 'PCL_Factor1');
modelPT_F1_PlotDV = modelPT_F1_Plot(1).YData';
modelPT_F1_PlotDV_Endpoints = modelPT_F1_Plot(2).YData';
modelPT_F1_ID_Overlap = strcmp(modelPT_F1.CoefficientNames, 'PCL_Factor1');

delete(modelPT_F1_Plot(1))
delete(modelPT_F1_Plot(2))
legend off

scatter(spreadRng - 10, modelPT_F1_PlotDV, 115, [0.5 0.5 0.5], 'filled', 'MarkerFaceAlpha', 0.5)
modelPT_F1_Scatter = scatter(data.PCL_Factor1 + spreadRng, modelPT_F1_PlotDV, 115, [0.5 0.5 0.5], 'filled', 'MarkerEdgeColor', 'none', 'LineWidth', 1, 'MarkerFaceAlpha', 0.65);
modelPT_F1_FromZero = sqrt(data.PCL_Factor1 .^ 2 + modelPT_F1_PlotDV .^ 2);
modelPT_F1_Scatter.AlphaData = modelPT_F1_FromZero;
modelPT_F1_Scatter.MarkerFaceAlpha = 'flat';

modelPT_F1_Plot = plotAdjustedResponse(modelPT_F1, 'PCL_Factor1');

delete(modelPT_F1_Plot(1))
modelPT_F1_Plot(2).LineStyle = '--';
modelPT_F1_Plot(2).LineWidth = 4;
modelPT_F1_Plot(2).Color = 'black';

set(gca, 'FontSize', 18)
title('IRI-PT ~ PCL-R F1', 'FontSize', 25)
subtitle({sprintf('β_Z = %.0s, P = %.3f, R^2_a_d_j = %.2f', modelPT_F1.Coefficients.Estimate(modelPT_F1_ID_Overlap), modelPT_F1.Coefficients.pValue(modelPT_F1_ID_Overlap), modelPT_F1.Rsquared.Adjusted)}, 'FontSize', 15)
xlabel('Interpersonal/Affective [Z]')
xlim([-3.25 3.25])
xticks([-3 0 3])
ylabel('Perspective Taking [Z]')
ylim([-3.25 3.25])
yticks([-3 0 3])
axis square
legend off

text(0.5, 0.95, sprintf('N = %.0f', sum(~isnan(data.PCL_Factor1))), ...
     'Units', 'normalized', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
     'FontSize', 15, 'BackgroundColor', 'white', 'EdgeColor', 'black');

% 2. IRI-EC ~ PCL-R F1

nexttile

line([0 0], [-3.5 3.5], 'LineWidth', 1.5, 'Color', [0.85 0.85 0.85], 'LineStyle', '--')
line([-3.5 3.5], [0 0], 'LineWidth', 1.5, 'Color', [0.85 0.85 0.85], 'LineStyle', '--')
hold on

rng default
spreadRng = 0.055 * randn(height(data), 1);

modelEC_F1 = fitlm(data, 'IRI_EC ~ PCL_Factor1 + Age + IQ', 'RobustOpts', 'on');
modelEC_F1_Plot = plotAdjustedResponse(modelEC_F1, 'PCL_Factor1');
modelEC_F1_PlotDV = modelEC_F1_Plot(1).YData';
modelEC_F1_PlotDV_Endpoints = modelEC_F1_Plot(2).YData';
modelEC_F1_ID_Overlap = strcmp(modelEC_F1.CoefficientNames, 'PCL_Factor1');

delete(modelEC_F1_Plot(1))
delete(modelEC_F1_Plot(2))
legend off

scatter(spreadRng - 10, modelEC_F1_PlotDV, 115, [0.5 0.5 0.5], 'filled', 'MarkerFaceAlpha', 0.5)
modelEC_F1_Scatter = scatter(data.PCL_Factor1 + spreadRng, modelEC_F1_PlotDV, 115, [0.5 0.5 0.5], 'filled', 'MarkerEdgeColor', 'none', 'LineWidth', 1, 'MarkerFaceAlpha', 0.65);
modelEC_F1_FromZero = sqrt(data.PCL_Factor1 .^ 2 + modelEC_F1_PlotDV .^ 2);
modelEC_F1_Scatter.AlphaData = modelEC_F1_FromZero;
modelEC_F1_Scatter.MarkerFaceAlpha = 'flat';

modelEC_F1_CI = fitlm(data.PCL_Factor1, modelEC_F1_PlotDV, 'RobustOpts', 'on');
modelEC_F1_CI_Plot = plot(modelEC_F1_CI, 'Marker', 'none');
delete(modelEC_F1_CI_Plot(1))
delete(modelEC_F1_CI_Plot(2))
modelEC_F1_CI_Plot(3).LineWidth = 2;
modelEC_F1_CI_Plot(4).LineWidth = 2;
modelEC_F1_CI_Plot(3).Color = 'green';
modelEC_F1_CI_Plot(4).Color = 'green';

modelEC_F1_Plot = plotAdjustedResponse(modelEC_F1, 'PCL_Factor1');

delete(modelEC_F1_Plot(1))
modelEC_F1_Plot(2).LineWidth = 4;
modelEC_F1_Plot(2).Color = 'black';

set(gca, 'FontSize', 18)
title('IRI-EC ~ PCL-R F1', 'FontSize', 25)
subtitle({sprintf('β_Z = %.2f, P_B_o_n = %.3f, R^2_a_d_j = %.2f', modelEC_F1.Coefficients.Estimate(modelEC_F1_ID_Overlap), modelEC_F1.Coefficients.pValue(modelEC_F1_ID_Overlap) * 2, modelEC_F1.Rsquared.Adjusted)}, 'FontSize', 15)
xlabel('Interpersonal/Affective [Z]')
xlim([-3.25 3.25])
xticks([-3 0 3])
ylabel('Empathic Concern [Z]')
ylim([-3.25 3.25])
yticks([-3 0 3])
axis square
legend off

text(0.5, 0.95, sprintf('N = %.0f', sum(~isnan(data.PCL_Factor1))), ...
     'Units', 'normalized', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
     'FontSize', 15, 'BackgroundColor', 'white', 'EdgeColor', 'black');

% 3. IRI-PT ~ PCL-R F2

nexttile

line([0 0], [-3.5 3.5], 'LineWidth', 1.5, 'Color', [0.85 0.85 0.85], 'LineStyle', '--')
line([-3.5 3.5], [0 0], 'LineWidth', 1.5, 'Color', [0.85 0.85 0.85], 'LineStyle', '--')
hold on

rng default
spreadRng = 0.055 * randn(height(dataF2), 1);

modelPT_F2 = fitlm(dataF2, 'IRI_PT ~ PCL_Factor2 + Age + IQ', 'RobustOpts', 'on');
modelPT_F2_Plot = plotAdjustedResponse(modelPT_F2, 'PCL_Factor2');
modelPT_F2_PlotDV = modelPT_F2_Plot(1).YData';
modelPT_F2_PlotDV_Endpoints = modelPT_F2_Plot(2).YData';
modelPT_F2_ID_Overlap = strcmp(modelPT_F2.CoefficientNames, 'PCL_Factor2');

delete(modelPT_F2_Plot(1))
delete(modelPT_F2_Plot(2))
legend off

scatter(spreadRng - 10, modelPT_F2_PlotDV, 115, [0.5 0.5 0.5], 'filled', 'MarkerFaceAlpha', 0.5)
modelPT_F2_Scatter = scatter(dataF2.PCL_Factor2 + spreadRng, modelPT_F2_PlotDV, 115, [0.5 0.5 0.5], 'filled', 'MarkerEdgeColor', 'none', 'LineWidth', 1, 'MarkerFaceAlpha', 0.65);
modelPT_F2_FromZero = sqrt(dataF2.PCL_Factor2 .^ 2 + modelPT_F2_PlotDV .^ 2);
modelPT_F2_Scatter.AlphaData = modelPT_F2_FromZero;
modelPT_F2_Scatter.MarkerFaceAlpha = 'flat';

modelPT_F2_CI = fitlm(dataF2.PCL_Factor2, modelPT_F2_PlotDV, 'RobustOpts', 'on');
modelPT_F2_CI_Plot = plot(modelPT_F2_CI, 'Marker', 'none');
delete(modelPT_F2_CI_Plot(1))
delete(modelPT_F2_CI_Plot(2))
modelPT_F2_CI_Plot(3).LineWidth = 2;
modelPT_F2_CI_Plot(4).LineWidth = 2;
modelPT_F2_CI_Plot(3).Color = 'green';
modelPT_F2_CI_Plot(4).Color = 'green';

modelPT_F2_Plot = plotAdjustedResponse(modelPT_F2, 'PCL_Factor2');

delete(modelPT_F2_Plot(1))
modelPT_F2_Plot(2).LineWidth = 4;
modelPT_F2_Plot(2).Color = 'black';

set(gca, 'FontSize', 18)
title('IRI-PT ~ PCL-R F2', 'FontSize', 25)
subtitle({sprintf('β_Z = %.2f, P_B_o_n = %.0s, R^2_a_d_j = %.2f', modelPT_F2.Coefficients.Estimate(modelPT_F2_ID_Overlap), modelPT_F2.Coefficients.pValue(modelPT_F2_ID_Overlap) * 2, modelPT_F2.Rsquared.Adjusted)}, 'FontSize', 15)
xlabel('Lifestyle/Antisocial [Z]')
xlim([-3.25 3.25])
xticks([-3 0 3])
ylabel('Perspective Taking [Z]')
ylim([-3.25 3.25])
yticks([-3 0 3])
axis square
legend off

text(0.5, 0.95, sprintf('N = %.0f', sum(~isnan(data.PCL_Factor2))), ...
     'Units', 'normalized', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
     'FontSize', 15, 'BackgroundColor', 'white', 'EdgeColor', 'black');

% 4. IRI-EC ~ PCL-R F2

nexttile

line([0 0], [-3.5 3.5], 'LineWidth', 1.5, 'Color', [0.85 0.85 0.85], 'LineStyle', '--')
line([-3.5 3.5], [0 0], 'LineWidth', 1.5, 'Color', [0.85 0.85 0.85], 'LineStyle', '--')
hold on

rng default
spreadRng = 0.055 * randn(height(dataF2), 1);

modelEC_F2 = fitlm(dataF2, 'IRI_EC ~ PCL_Factor2 + Age + IQ', 'RobustOpts', 'on');
modelEC_F2_Plot = plotAdjustedResponse(modelEC_F2, 'PCL_Factor2');
modelEC_F2_PlotDV = modelEC_F2_Plot(1).YData';
modelEC_F2_PlotDV_Endpoints = modelEC_F2_Plot(2).YData';
modelEC_F2_ID_Overlap = strcmp(modelEC_F2.CoefficientNames, 'PCL_Factor2');

delete(modelEC_F2_Plot(1))
delete(modelEC_F2_Plot(2))
legend off

scatter(spreadRng - 10, modelEC_F2_PlotDV, 115, [0.5 0.5 0.5], 'filled', 'MarkerFaceAlpha', 0.5)
modelEC_F2_Scatter = scatter(dataF2.PCL_Factor2 + spreadRng, modelEC_F2_PlotDV, 115, [0.5 0.5 0.5], 'filled', 'MarkerEdgeColor', 'none', 'LineWidth', 1, 'MarkerFaceAlpha', 0.65);
modelEC_F2_FromZero = sqrt(dataF2.PCL_Factor2 .^ 2 + modelEC_F2_PlotDV .^ 2);
modelEC_F2_Scatter.AlphaData = modelEC_F2_FromZero;
modelEC_F2_Scatter.MarkerFaceAlpha = 'flat';

modelEC_F2_CI = fitlm(dataF2.PCL_Factor2, modelEC_F2_PlotDV, 'RobustOpts', 'on');
modelEC_F2_CI_Plot = plot(modelEC_F2_CI, 'Marker', 'none');
delete(modelEC_F2_CI_Plot(1))
delete(modelEC_F2_CI_Plot(2))
modelEC_F2_CI_Plot(3).LineWidth = 2;
modelEC_F2_CI_Plot(4).LineWidth = 2;
modelEC_F2_CI_Plot(3).Color = 'green';
modelEC_F2_CI_Plot(4).Color = 'green';

modelEC_F2_Plot = plotAdjustedResponse(modelEC_F2, 'PCL_Factor2');

delete(modelEC_F2_Plot(1))
modelEC_F2_Plot(2).LineWidth = 4;
modelEC_F2_Plot(2).Color = 'black';

set(gca, 'FontSize', 18)
title('IRI-EC ~ PCL-R F2', 'FontSize', 25)
subtitle({sprintf('β_Z = %.2f, P_B_o_n = %.0s, R^2_a_d_j = %.2f', modelEC_F2.Coefficients.Estimate(modelEC_F2_ID_Overlap), modelEC_F2.Coefficients.pValue(modelEC_F2_ID_Overlap) * 2, modelEC_F2.Rsquared.Adjusted)}, 'FontSize', 15)
xlabel('Lifestyle/Antisocial [Z]')
xlim([-3.25 3.25])
xticks([-3 0 3])
ylabel('Empathic Concern [Z]')
ylim([-3.25 3.25])
yticks([-3 0 3])
axis square
legend off

text(0.5, 0.95, sprintf('N = %.0f', sum(~isnan(data.PCL_Factor2))), ...
     'Units', 'normalized', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
     'FontSize', 15, 'BackgroundColor', 'white', 'EdgeColor', 'black');

%% 1.3: IRI subscales ~ PCL-R total & PCL-R group

clc; clear; close all

load data

figure('units', 'normalized', 'outerposition', [0 0 1 1])
T = tiledlayout(2, 4);

% 1. PCL-R total

nexttile
H_PCL_Total = histogram(data.PCL_Total, 'BinWidth', 1, 'BinLimits', [-0.5 40.5], 'FaceColor', 'blue', 'EdgeColor', 'none', 'FaceAlpha', 1);
H_PCL_Total.BinLimits = [-0.5 40.5];
xlim([-1 41])
hold on

for i = 1 : length(H_PCL_Total.BinCounts)
    
     binHeight = H_PCL_Total.BinCounts(i);
     binCentre = (H_PCL_Total.BinEdges(i) + H_PCL_Total.BinEdges(i + 1)) / 2; 

     if binCentre <= 20
        barColor = [0.5 0.5 0.5];
     elseif binCentre >= 30
        barColor = 'red';
     else
        barColor = [0.85 0.85 0.85];
     end

     bar(binCentre, binHeight, 'FaceColor', barColor, 'EdgeColor', 'none', 'LineWidth', 1.5, 'BarWidth', 1);
     
end

set(gca, 'FontSize', 18)
title('PCL-R', 'FontSize', 25)
subtitle(sprintf('X^X_X N = %.0f X^X_X', sum(~isnan(data.PCL_Total))), 'FontSize', 15)
xlabel('Highpathy Checklist-Revised')
xlim([-1 41])
xticks([0 ...
     median(data.PCL_Total(data.Status == 'Low')) ...
     median(data.PCL_Total(data.Status == 'High')) ...
     40])
ylabel('Frequency')
ylim([0 max(H_PCL_Total.BinCounts)])
yticks([0 max(H_PCL_Total.BinCounts)])
box off
axis square

line([median(data.PCL_Total(data.Status == 'Low')) median(data.PCL_Total(data.Status == 'Low'))], [0 100], 'LineWidth', 2, 'Color', 'black', 'LineStyle', '--')
line([median(data.PCL_Total(data.Status == 'High')) median(data.PCL_Total(data.Status == 'High'))], [0 100], 'LineWidth', 2, 'Color', 'black', 'LineStyle', '--')

[F_PCL_Total, XI_PCL_Total] = ksdensity(data.PCL_Total, 'BandWidth', 2);
plot(XI_PCL_Total, F_PCL_Total * sum(H_PCL_Total.BinWidth * H_PCL_Total.BinCounts), 'Color', [0.5 0.5 0.5], 'LineWidth', 4);

legendInfo = get(gca, 'Children');
legend([legendInfo(end - 2) legendInfo(end - 40) legendInfo(2)], {'Low' 'High' 'Med'}, 'FontSize', 15, 'Location', 'NorthWest', 'box', 'off')

% 2. IRI-PT ~ Group

data.IRI_PT = zscore(data.IRI_PT);
data.IRI_EC = zscore(data.IRI_EC);

dataStatus.IRI_PT = zscore(dataStatus.IRI_PT);
dataStatus.IRI_EC = zscore(dataStatus.IRI_EC);

nexttile

line([0 3.5], [0 0], 'LineWidth', 1.5, 'Color', [0.85 0.85 0.85], 'LineStyle', '--')
hold on

rng default
spreadLow = 0.055 * randn(sum(dataStatus.Status == 'Low'), 1);
spreadHigh = 0.055 * randn(sum(dataStatus.Status == 'High'), 1);

modelPT = fitlm(dataStatus, 'IRI_PT ~ Status + Age + IQ', 'RobustOpts', 'on');
modelPT.coefCI;
modelPT_Plot = plotAdjustedResponse(modelPT, 'Status');
modelPT_PlotDV = modelPT_Plot(1).YData';
modelPT_PlotDV_Endpoints = modelPT_Plot(2).YData';
modelPT_ID = strcmp(modelPT.CoefficientNames, 'Status_High');

delete(modelPT_Plot(1))
delete(modelPT_Plot(2))
legend off

colours = [0.5 0.5 0.5; [228 26 28] / 255];
VP_1 = violinplot(modelPT_PlotDV, dataStatus.Status, ...
     'ShowMean', false, ...
     'ShowMedian', false, ...
     'ShowData', false, ...
     'ShowBox', false, ...
     'ViolinColor', colours, ...
     'BoxColor', [1 1 1], ...
     'HalfViolin', 'right', ...
     'ViolinAlpha', 0, ...
     'Width', 0.385);
for i = 1 : 2
     VP_1(i).BoxColor = colours(i, :);
     VP_1(i).EdgeColor = [1 1 1];
     VP_1(i).EdgeColor = colours(i, :);
     VP_1(i).ViolinPlot.LineWidth = 2;
end

scatter(spreadLow + 1, modelPT_PlotDV(dataStatus.Status == 'Low'), 115, [0.5 0.5 0.5], 'filled', 'MarkerFaceAlpha', 0.5);
scatter(spreadHigh + 2, modelPT_PlotDV(dataStatus.Status == 'High'), 115, [228 26 28] / 255, 'filled', 'MarkerFaceAlpha', 0.5);

modelPT_CI = fitlm(dataStatus.Status, modelPT_PlotDV, 'RobustOpts', 'on');
modelPT_CI_Plot = plot(modelPT_CI, 'Marker', 'none');
delete(modelPT_CI_Plot(1))
delete(modelPT_CI_Plot(2))
modelPT_CI_Plot(3).LineWidth = 2;
modelPT_CI_Plot(4).LineWidth = 2;
modelPT_CI_Plot(3).Color = 'green';
modelPT_CI_Plot(4).Color = 'green';

modelPT_Plot = plotAdjustedResponse(modelPT, 'Status');

delete(modelPT_Plot(1))
modelPT_Plot(2).LineWidth = 4;
modelPT_Plot(2).Color = 'black';

modelPT_ErrorBarLow = errorbar(1, modelPT_PlotDV_Endpoints(1), std(modelPT_PlotDV(dataStatus.Status == 'Low')));
modelPT_ErrorBarLow.LineWidth = 4;
modelPT_ErrorBarLow.Color = 'black';

modelPT_ErrorBarHigh = errorbar(2, modelPT_PlotDV_Endpoints(2), std(modelPT_PlotDV(dataStatus.Status == 'High')));
modelPT_ErrorBarHigh.LineWidth = 4;
modelPT_ErrorBarHigh.Color = 'black';

set(gca, 'FontSize', 18)
title('IRI-PT ~ Group', 'FontSize', 25)
subtitle({sprintf('P_B_o_n = %.3f, R^2_a_d_j = %.2f', modelPT.Coefficients.pValue(modelPT_ID) * 2, modelPT.Rsquared.Adjusted)}, 'FontSize', 15)
xlabel('Highpathy group')
xlim([0.5 2.5])
xticks([1 2])
xticklabels({'Low' 'High'})
ylabel('Perspective Taking [Z]')
ylim([-3 3])
yticks([-3 0 3])
box off
axis square
legend off

modelPT = fitlm(dataStatus, 'IRI_PT ~ Status + Age + IQ', 'RobustOpts', 'on');
modelPT_ID = strcmp(modelPT.CoefficientNames, 'Status_High');
modelPT_Beta = modelPT.Coefficients.Estimate(modelPT_ID);

modelPT_Residuals = fitlm(dataStatus, 'IRI_PT ~ Age + IQ', 'RobustOpts', 'on');
modelPT_CohensD = computeCohen_d(modelPT_Residuals.Residuals.Raw(dataStatus.Status == 'High'), modelPT_Residuals.Residuals.Raw(dataStatus.Status == 'Low'));

modelPT_Text = "Cohen's D =";
modelPT_Text2 = sprintf('%.2f', modelPT_CohensD);

text(0.5, 0.928, sprintf('%s', modelPT_Text), ...
     'Units', 'normalized', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
     'FontSize', 15);

text(0.5, 0.818, sprintf('%s', modelPT_Text2), ...
     'Units', 'normalized', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
     'FontSize', 15);

text(0.25, 0.05, sprintf('N = %.0f', sum(dataStatus.Status == 'Low')), ...
     'Units', 'normalized', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
     'Background', 'white', 'EdgeColor', 'black', ...
     'FontSize', 15);

text(0.75, 0.05, sprintf('N = %.0f', sum(dataStatus.Status == 'High')), ...
     'Units', 'normalized', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
     'Background', 'white', 'EdgeColor', 'black', ...
     'FontSize', 15);

% 3. IRI-EC ~ Status

nexttile

line([0 3.5], [0 0], 'LineWidth', 1.5, 'Color', [0.85 0.85 0.85], 'LineStyle', '--')
hold on

rng default
spreadLow = 0.055 * randn(sum(dataStatus.Status == 'Low'), 1);
spreadHigh = 0.055 * randn(sum(dataStatus.Status == 'High'), 1);

modelEC = fitlm(dataStatus, 'IRI_EC ~ Status + Age + IQ', 'RobustOpts', 'on');
modelEC.coefCI;
modelEC_Plot = plotAdjustedResponse(modelEC, 'Status');
modelEC_PlotDV = modelEC_Plot(1).YData';
modelEC_PlotDV_Endpoints = modelEC_Plot(2).YData';
modelEC_ID = strcmp(modelEC.CoefficientNames, 'Status_High');

delete(modelEC_Plot(1))
delete(modelEC_Plot(2))
legend off

colours = [0.5 0.5 0.5; [228 26 28] / 255];
VP_1 = violinplot(modelEC_PlotDV, dataStatus.Status, ...
     'ShowMean', false, ...
     'ShowMedian', false, ...
     'ShowData', false, ...
     'ShowBox', false, ...
     'ViolinColor', colours, ...
     'BoxColor', [1 1 1], ...
     'HalfViolin', 'right', ...
     'ViolinAlpha', 0, ...
     'Width', 0.385);
for i = 1 : 2
     VP_1(i).BoxColor = colours(i, :);
     VP_1(i).EdgeColor = [1 1 1];
     VP_1(i).EdgeColor = colours(i, :);
     VP_1(i).ViolinPlot.LineWidth = 2;
end

scatter(spreadLow + 1, modelEC_PlotDV(dataStatus.Status == 'Low'), 115, [0.5 0.5 0.5], 'filled', 'MarkerFaceAlpha', 0.5);
scatter(spreadHigh + 2, modelEC_PlotDV(dataStatus.Status == 'High'), 115, [228 26 28] / 255, 'filled', 'MarkerFaceAlpha', 0.5);

modelEC_CI = fitlm(dataStatus.Status, modelEC_PlotDV, 'RobustOpts', 'on');
modelEC_CI_Plot = plot(modelEC_CI, 'Marker', 'none');
delete(modelEC_CI_Plot(1))
delete(modelEC_CI_Plot(2))
modelEC_CI_Plot(3).LineWidth = 2;
modelEC_CI_Plot(4).LineWidth = 2;
modelEC_CI_Plot(3).Color = 'green';
modelEC_CI_Plot(4).Color = 'green';

modelEC_Plot = plotAdjustedResponse(modelEC, 'Status');

delete(modelEC_Plot(1))
modelEC_Plot(2).LineWidth = 4;
modelEC_Plot(2).Color = 'black';

modelEC_ErrorBarLow = errorbar(1, modelEC_PlotDV_Endpoints(1), std(modelEC_PlotDV(dataStatus.Status == 'Low')));
modelEC_ErrorBarLow.LineWidth = 4;
modelEC_ErrorBarLow.Color = 'black';

modelEC_ErrorBarHigh = errorbar(2, modelEC_PlotDV_Endpoints(2), std(modelEC_PlotDV(dataStatus.Status == 'High')));
modelEC_ErrorBarHigh.LineWidth = 4;
modelEC_ErrorBarHigh.Color = 'black';

set(gca, 'FontSize', 18)
title('IRI-EC ~ Group', 'FontSize', 25)
subtitle({sprintf('P_B_o_n = %.0s, R^2_a_d_j = %.2f', modelEC.Coefficients.pValue(modelEC_ID) * 2, modelEC.Rsquared.Adjusted)}, 'FontSize', 15)
xlabel('Highpathy group')
xlim([0.5 2.5])
xticks([1 2])
xticklabels({'Low' 'High'})
ylabel('Empathic Concern [Z]')
ylim([-3 3])
yticks([-3 0 3])
box off
axis square
legend off

modelEC = fitlm(dataStatus, 'IRI_EC ~ Status + Age + IQ', 'RobustOpts', 'on');
modelEC_ID = strcmp(modelEC.CoefficientNames, 'Status_High');
modelEC_Beta = modelEC.Coefficients.Estimate(modelEC_ID);

modelEC_Residuals = fitlm(dataStatus, 'IRI_EC ~ Age + IQ', 'RobustOpts', 'on');
modelEC_CohensD = computeCohen_d(modelEC_Residuals.Residuals.Raw(dataStatus.Status == 'High'), modelEC_Residuals.Residuals.Raw(dataStatus.Status == 'Low'));

modelEC_Text = "Cohen's D =";
modelEC_Text2 = sprintf('%.2f', modelEC_CohensD);

text(0.5, 0.928, sprintf('%s', modelEC_Text), ...
     'Units', 'normalized', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
     'FontSize', 15);

text(0.5, 0.818, sprintf('%s', modelEC_Text2), ...
     'Units', 'normalized', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
     'FontSize', 15);

text(0.25, 0.05, sprintf('N = %.0f', sum(dataStatus.Status == 'Low')), ...
     'Units', 'normalized', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
     'Background', 'white', 'EdgeColor', 'black', ...
     'FontSize', 15);

text(0.75, 0.05, sprintf('N = %.0f', sum(dataStatus.Status == 'High')), ...
     'Units', 'normalized', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
     'Background', 'white', 'EdgeColor', 'black', ...
     'FontSize', 15);

% 4. Sensitivity: IRI ~ Group

modelPT_Residuals = fitlm(dataStatus, 'IRI_PT ~ Age + IQ', 'RobustOpts', 'on');
modelPT_CohensD = computeCohen_d(modelPT_Residuals.Residuals.Raw(dataStatus.Status == 'High'), modelPT_Residuals.Residuals.Raw(dataStatus.Status == 'Low'));

modelEC_Residuals = fitlm(dataStatus, 'IRI_EC ~ Age + IQ', 'RobustOpts', 'on');
modelEC_CohensD = computeCohen_d(modelEC_Residuals.Residuals.Raw(dataStatus.Status == 'High'), modelEC_Residuals.Residuals.Raw(dataStatus.Status == 'Low'));

modelPT_EC_Residuals = fitlm(dataStatus, 'IRI_PT ~ Age + IQ + IRI_EC', 'RobustOpts', 'on');
modelPT_EC_CohensD = computeCohen_d(modelPT_EC_Residuals.Residuals.Raw(dataStatus.Status == 'High'), modelPT_EC_Residuals.Residuals.Raw(dataStatus.Status == 'Low'));

modelEC_PT_Residuals = fitlm(dataStatus, 'IRI_EC ~ Age + IQ + IRI_PT', 'RobustOpts', 'on');
modelEC_PT_CohensD = computeCohen_d(modelEC_PT_Residuals.Residuals.Raw(dataStatus.Status == 'High'), modelEC_PT_Residuals.Residuals.Raw(dataStatus.Status == 'Low'));

CohensDs = [modelPT_CohensD modelEC_CohensD modelPT_EC_CohensD modelEC_PT_CohensD];

nexttile

bar(1, CohensDs(1), 'LineWidth', 2, 'EdgeColor', 'black', 'FaceColor', 'w', 'BarWidth', 0.5)
hold on
bar(2, CohensDs(2), 'LineWidth', 2, 'EdgeColor', 'black', 'FaceColor', 'w', 'BarWidth', 0.5)
hold on
bar(3, CohensDs(3), 'LineWidth', 2, 'EdgeColor', 'black', 'FaceColor', 'w', 'BarWidth', 0.5)
hold on
bar(4, CohensDs(4), 'LineWidth', 2, 'EdgeColor', 'black', 'FaceColor', 'w', 'BarWidth', 0.5)

set(gca, 'FontSize', 18)
title('IRI ~ Group', 'FontSize', 25)
subtitle(sprintf('X^X_X N = %.0f vs %.0f X^X_X', sum(dataStatus.Status == 'High'), sum(dataStatus.Status == 'Low')), 'FontSize', 15)
xlim([0 5])
xticks([1 : 4])
xticklabels({'IRI-PT' 'IRI-EC' 'PT | EC' 'EC | PT'})
xtickangle(25)
ylabel("Cohen's D: High vs low")
perc = max(abs(CohensDs))*0.15;
ylim([min(CohensDs)-perc 0])
yticks([min(CohensDs) 0])
yticklabels({sprintf('%.2f', min(CohensDs)) 0})
axis square
box off

text(1/5, 0.05, '*', ...
     'Units', 'normalized', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
     'Color', 'black', ...
     'FontSize', 25);

text(2/5, 0.05, '*', ...
     'Units', 'normalized', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
     'Color', 'black', ...
     'FontSize', 25);

text(4/5, 0.05, '*', ...
     'Units', 'normalized', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
     'Color', 'black', ...
     'FontSize', 25);

scatter(-10, -10)
tileInfo = get(gca, 'Children');
legend(tileInfo(1), 'P_B_o_n < 0.05', 'Location', 'layout', 'Orientation', 'horizontal', 'FontSize', 15);

%% 1.4 Correlation matrices 1

clc; clear; close all

load data

dataStatus.Status_Num = dataStatus.Status;
dataStatus.Status_Num(dataStatus.Status_Num == 'High') = '1';
dataStatus.Status_Num(dataStatus.Status_Num ~= '1') = '0';
dataStatus.Status_Num = string(dataStatus.Status_Num);
dataStatus.Status_Num = str2double(dataStatus.Status_Num);

% 1. Total

arrayTotal = table2array(data(:, {'IRI_PT' 'IRI_EC' 'PCL_Total' 'PCL_Factor1' 'PCL_Factor2' 'Age' 'IQ' 'eTIV'}));
namesTotal = {'IRI-PT' 'IRI-EC' 'PCL-R' 'PCL-R F1' 'PCL-R F2' 'Age' 'IQ' 'TIV'};
tableTotal = array2table(arrayTotal, 'VariableNames', namesTotal);

arrayTotal = table2array(tableTotal);
[rhosTotal, P_ValuesTotal] = corr(arrayTotal, 'rows', 'pairwise', 'type', 'Spearman');
rhosTotal(logical(eye(size(rhosTotal)))) = 0;

matrixTotal = triu(rhosTotal) + tril(P_ValuesTotal);
matrixTotal = matrixTotal .* ~eye(size(matrixTotal));
matrixTotal = array2table(matrixTotal, 'VariableNames', namesTotal);
matrixTotal.Row = namesTotal

% 2. Low

arrayLow = table2array(dataLow(:, {'IRI_PT' 'IRI_EC' 'PCL_Total' 'PCL_Factor1' 'PCL_Factor2' 'Age' 'IQ' 'eTIV'}));
namesLow = {'IRI-PT' 'IRI-EC' 'PCL-R' 'PCL-R F1' 'PCL-R F2' 'Age' 'IQ' 'TIV'};
tableLow = array2table(arrayLow, 'VariableNames', namesLow);

arrayLow = table2array(tableLow);
[rhosLow, P_ValuesLow] = corr(arrayLow, 'rows', 'pairwise', 'type', 'Spearman');
rhosLow(logical(eye(size(rhosLow)))) = 0;

matrixLow = triu(rhosLow) + tril(P_ValuesLow);
matrixLow = matrixLow .* ~eye(size(matrixLow));
matrixLow = array2table(matrixLow, 'VariableNames', namesLow);
matrixLow.Row = namesLow

% 3. High

arrayHigh = table2array(dataHigh(:, {'IRI_PT' 'IRI_EC' 'PCL_Total' 'PCL_Factor1' 'PCL_Factor2' 'Age' 'IQ' 'eTIV'}));
namesHigh = {'IRI-PT' 'IRI-EC' 'PCL-R' 'PCL-R F1' 'PCL-R F2' 'Age' 'IQ' 'TIV'};
tableHigh = array2table(arrayHigh, 'VariableNames', namesHigh);

arrayHigh = table2array(tableHigh);
[rhosHigh, P_ValuesHigh] = corr(arrayHigh, 'rows', 'pairwise', 'type', 'Spearman');
rhosHigh(logical(eye(size(rhosHigh)))) = 0;

matrixHigh = triu(rhosHigh) + tril(P_ValuesHigh);
matrixHigh = matrixHigh .* ~eye(size(matrixHigh));
matrixHigh = array2table(matrixHigh, 'VariableNames', namesHigh);
matrixHigh.Row = namesHigh

save correlationMatrices

%% 1.5 Correlation matrices 2

clc; clear; close all

load correlationMatrices

N_Tests = height(P_ValuesTotal) * (height(P_ValuesTotal) - 1) / 2;

P_ValuesTotalBon = P_ValuesTotal * N_Tests;
P_ValuesTotalBon(P_ValuesTotalBon > 1) = 1;

P_ValuesLowBon = P_ValuesLow * N_Tests;
P_ValuesLowBon(P_ValuesLowBon > 1) = 1;

P_ValuesHighBon = P_ValuesHigh * N_Tests;
P_ValuesHighBon(P_ValuesHighBon > 1) = 1;

figure('units', 'normalized', 'outerposition', [0 0 1 1])
T = tiledlayout(5, 3);

% 1. Total

nexttile([2 1])
imagesc(rhosTotal, [-0.5 0.5]);
colormap(RdBu_r)
xticks([1 : width(namesTotal)])
xticklabels(namesTotal)
xtickangle(45)
yticks([1 : width(namesTotal)])
yticklabels(namesTotal)
set(gca, 'FontSize', 18)
title('Total', 'FontSize', 25)
subtitle(sprintf('N = %.0f', height(data)), 'FontSize', 20)
box off

for i = 1 : width(namesTotal)
     
     for k = 1 : width(namesTotal)

     if logical(P_ValuesTotalBon(i, k) < 0.05) & abs(rhosTotal(i, k)) >= 0.3
     text(i, k, sprintf('%.2f', rhosTotal(i, k)), 'Color', 'w', 'FontSize', 13, 'HorizontalAlignment', 'center');
     else
     end

     if logical(P_ValuesTotalBon(i, k) < 0.05) & abs(rhosTotal(i, k)) < 0.3
     text(i, k, sprintf('%.2f', rhosTotal(i, k)), 'Color', 'black', 'FontSize', 13, 'HorizontalAlignment', 'center');
     else
     end

     end
     
end

% 2. Low

nexttile([2 1])
imagesc(rhosLow, [-0.5 0.5]);
colormap(RdBu_r)
xticks([1 : width(namesLow)])
xticklabels(namesLow)
xtickangle(45)
yticks([1 : width(namesLow)])
yticklabels(namesLow)
set(gca, 'FontSize', 18)
title('Low', 'FontSize', 25)
subtitle(sprintf('N = %.0f', height(dataLow)), 'FontSize', 20)
box off

for i = 1 : width(namesLow)
     
     for k = 1 : width(namesLow)

     if logical(P_ValuesLowBon(i, k) < 0.05) & abs(rhosLow(i, k)) >= 0.3
     text(i, k, sprintf('%.2f', rhosLow(i, k)), 'Color', 'w', 'FontSize', 13, 'HorizontalAlignment', 'center');
     else
     end

     if logical(P_ValuesLowBon(i, k) < 0.05) & abs(rhosLow(i, k)) < 0.3
     text(i, k, sprintf('%.2f', rhosLow(i, k)), 'Color', 'black', 'FontSize', 13, 'HorizontalAlignment', 'center');
     else
     end

     end
     
end

% 3. High

nexttile([2 1])
imagesc(rhosHigh, [-0.5 0.5]);
colormap(RdBu_r)
xticks([1 : width(namesHigh)])
xticklabels(namesHigh)
xtickangle(45)
yticks([1 : width(namesHigh)])
yticklabels(namesHigh)
set(gca, 'FontSize', 18)
title('High', 'FontSize', 25)
subtitle(sprintf('N = %.0f', height(dataHigh)), 'FontSize', 20)
box off

for i = 1 : width(namesHigh)
     
     for k = 1 : width(namesHigh)

     if logical(P_ValuesHighBon(i, k) < 0.05) & abs(rhosHigh(i, k)) >= 0.3
     text(i, k, sprintf('%.2f', rhosHigh(i, k)), 'Color', 'w', 'FontSize', 13, 'HorizontalAlignment', 'center');
     else
     end

     if logical(P_ValuesHighBon(i, k) < 0.05) & abs(rhosHigh(i, k)) < 0.3
     text(i, k, sprintf('%.2f', rhosHigh(i, k)), 'Color', 'black', 'FontSize', 13, 'HorizontalAlignment', 'center');
     else
     end

     end
     
end

bar = colorbar('FontSize', 20);
bar.Label.String = "Spearman's ρ";
bar.Ticks = [-0.5 0 0.5];