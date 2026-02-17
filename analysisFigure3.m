%% Figure 3

%% 1.1 Prepare the training and test sets: IRI-PT, IRI-EC, and PCL-R F1 
% (and PCL-R total score, for the supplement)

clc; clear; close all

load data

% Select variables of interest
dataNonSA = data(:, {'URSI' 'IRI_PT' 'IRI_EC' 'PCL_Total' 'PCL_Factor1' 'Age' 'IQ' 'eTIV'});
dataSA = data(:, labelsGlasserSA);
data = [dataNonSA dataSA];

% Partition X data into the training and test sets
rng default
CV = cvpartition(size(data, 1), 'HoldOut', 0.2);
dataX_Train = data(CV.training, :);
dataX_Test = data(CV.test, :);

% Correct X data for covariates in the training set, to prevent data
% leakage
measure = 'Area';
covariates = '~ Age + IQ + eTIV';
measureOfInterest = strfind(dataX_Train.Properties.VariableNames, measure);
measureOfInterestID = find(~cellfun(@isempty, measureOfInterest));

for parcelID = 1 : numel(measureOfInterestID)

    parcel = dataX_Train.Properties.VariableNames{measureOfInterestID(parcelID)};
    fprintf('In the training set, correcting: %s\n', parcel)
    modelToBeFitted = sprintf('%s %s', parcel, covariates);
    model = fitlm(dataX_Train, modelToBeFitted, 'RobustOpts', 'on');
    X_Train(:, parcelID) = model.Residuals.Raw;

end
dataX_Train(:, labelsGlasserSA) = array2table(X_Train);

% Correct X data for covariates in the test set
measure = 'Area';
covariates = '~ Age + IQ + eTIV';
measureOfInterest = strfind(dataX_Test.Properties.VariableNames, measure);
measureOfInterestID = find(~cellfun(@isempty, measureOfInterest));

for parcelID = 1 : numel(measureOfInterestID)

    parcel = dataX_Test.Properties.VariableNames{measureOfInterestID(parcelID)};
    fprintf('In the test set, correcting: %s\n', parcel)
    modelToBeFitted = sprintf('%s %s', parcel, covariates);
    model = fitlm(dataX_Test, modelToBeFitted, 'RobustOpts', 'on');
    X_Test(:, parcelID) = model.Residuals.Raw;

end
dataX_Test(:, labelsGlasserSA) = array2table(X_Test);

clearvars -except CV X* Y* data dataX* dataY*

% save dataTrainTestSA

%% 1.2 Prepare the training and test sets: PCL-R F2

clc; clear; close all

load data

% Select variables of interest
dataNonSA = data(:, {'URSI' 'PCL_Factor2' 'Age' 'IQ' 'eTIV'});
dataSA = data(:, labelsGlasserSA);
data = [dataNonSA dataSA];

% Specify Y
Y = data.PCL_Factor2;

% Drop missing observations
[Y, IDX] = rmmissing(Y);
data(IDX, :) = [];

% Partition X data into the training and test sets
rng default
CV = cvpartition(size(data, 1), 'HoldOut', 0.2);
dataX_Train = data(CV.training, :);
dataX_Test = data(CV.test, :);

% Correct X data for covariates in the training set, to prevent data
% leakage
measure = 'Area';
covariates = '~ Age + IQ + eTIV';
measureOfInterest = strfind(dataX_Train.Properties.VariableNames, measure);
measureOfInterestID = find(~cellfun(@isempty, measureOfInterest));

for parcelID = 1 : numel(measureOfInterestID)

    parcel = dataX_Train.Properties.VariableNames{measureOfInterestID(parcelID)};
    fprintf('In the training set, correcting: %s\n', parcel)
    modelToBeFitted = sprintf('%s %s', parcel, covariates);
    model = fitlm(dataX_Train, modelToBeFitted, 'RobustOpts', 'on');
    X_Train(:, parcelID) = model.Residuals.Raw;

end
dataX_Train(:, labelsGlasserSA) = array2table(X_Train);

% Correct X data for covariates in the test set
measure = 'Area';
covariates = '~ Age + IQ + eTIV';
measureOfInterest = strfind(dataX_Test.Properties.VariableNames, measure);
measureOfInterestID = find(~cellfun(@isempty, measureOfInterest));

for parcelID = 1 : numel(measureOfInterestID)

    parcel = dataX_Test.Properties.VariableNames{measureOfInterestID(parcelID)};
    fprintf('In the test set, correcting: %s\n', parcel)
    modelToBeFitted = sprintf('%s %s', parcel, covariates);
    model = fitlm(dataX_Test, modelToBeFitted, 'RobustOpts', 'on');
    X_Test(:, parcelID) = model.Residuals.Raw;

end
dataX_Test(:, labelsGlasserSA) = array2table(X_Test);

clearvars -except CV X* data dataX* dataY*

% save dataTrainTestPCL_Factor2_SA

%% 2.1 Fit the model in the training set: IRI-PT

clc; clear; close all

load dataTrainTestSA

% Specify X and Y data
X_TrainIRI_PT = X_Train;
X_TestIRI_PT = X_Test;
Y_IRI_PT = data.IRI_PT;

% Partition Y data
Y_TrainIRI_PT = Y_IRI_PT(CV.training, :);
Y_TestIRI_PT = Y_IRI_PT(CV.test, :);

clearvars -except *IRI_PT

% Normalise X data
a = 0;
b = 1;

minX_TrainIRI_PT = min(X_TrainIRI_PT);
maxX_TrainIRI_PT = max(X_TrainIRI_PT);

X_TrainIRI_PT = (b - a) .* ((X_TrainIRI_PT - minX_TrainIRI_PT) ./ (maxX_TrainIRI_PT - minX_TrainIRI_PT)) + a;
X_TestIRI_PT = (b - a) .* ((X_TestIRI_PT - minX_TrainIRI_PT) ./ (maxX_TrainIRI_PT - minX_TrainIRI_PT)) + a;

% Specify N for CV folds
N_Folds = 10;

% Specify 1,000 lambdas
lambdas = logspace(-3, 0, 1000);

% Fit ridge regression
X_TrainTrans = X_TrainIRI_PT';
rng default
modelIRI_PT = fitrlinear(X_TrainTrans, Y_TrainIRI_PT, 'ObservationsIn', 'columns', 'KFold', N_Folds, 'Lambda', lambdas, 'Learner', 'leastsquares', 'Regularization', 'ridge', 'Solver', 'lbfgs');
MSEsIRI_PT = kfoldLoss(modelIRI_PT);

% Fit the same ridge regression but without CV
modelNonCV_IRI_PT = fitrlinear(X_TrainTrans, Y_TrainIRI_PT, 'ObservationsIn', 'columns', 'Lambda', lambdas, 'Learner', 'leastsquares', 'Regularization', 'ridge', 'Solver', 'lbfgs');

% Extract the non-CV betas
modelNonCV_BetasIRI_PT = modelNonCV_IRI_PT.Beta;

% Select the modelIRI_PT corresponding to the minimal MSE
[MSE_MinIRI_PT, MSE_MinID_IRI_PT] = min(MSEsIRI_PT);
modelFinalIRI_PT = selectModels(modelNonCV_IRI_PT, MSE_MinID_IRI_PT);

% Extract the final betas and intercept
modelFinalBetasIRI_PT = modelFinalIRI_PT.Beta;
modelFinalInterceptIRI_PT = modelFinalIRI_PT.Bias;

% Predict Y in the test set
Y_TestHatIRI_PT = cat(2, ones(size(X_TestIRI_PT, 1), 1), X_TestIRI_PT) * [modelFinalInterceptIRI_PT; modelFinalBetasIRI_PT];

% Calculate stats
L2_NormIRI_PT = norm(modelFinalBetasIRI_PT); % L2 norm
SSE_IRI_PT = sum((Y_TestIRI_PT - Y_TestHatIRI_PT) .^ 2); % Sum of Squared Errors
SST_IRI_PT = sum((Y_TestIRI_PT - mean(Y_TestIRI_PT)) .^ 2); % Total Sum of Squares
MSE_IRI_PT = mean((Y_TestIRI_PT - Y_TestHatIRI_PT) .^ 2); % Mean Squared Error
MAE_IRI_PT = mean(abs(Y_TestIRI_PT - Y_TestHatIRI_PT)); % Mean Absolute Error
R_SquaredIRI_PT = 1 - (SSE_IRI_PT / SST_IRI_PT) % Coefficient of Determination

% save ridgeIRI_PT_SA

%% 2.2 Fit the model in the training set: IRI-EC

clc; clear; close all

load dataTrainTestSA

% Specify X and Y data
X_TrainIRI_EC = X_Train;
X_TestIRI_EC = X_Test;
Y_IRI_EC = data.IRI_EC;

% Partition Y data
Y_TrainIRI_EC = Y_IRI_EC(CV.training, :);
Y_TestIRI_EC = Y_IRI_EC(CV.test, :);

clearvars -except *IRI_EC

% Normalise X data
a = 0;
b = 1;

minX_TrainIRI_EC = min(X_TrainIRI_EC);
maxX_TrainIRI_EC = max(X_TrainIRI_EC);

X_TrainIRI_EC = (b - a) .* ((X_TrainIRI_EC - minX_TrainIRI_EC) ./ (maxX_TrainIRI_EC - minX_TrainIRI_EC)) + a;
X_TestIRI_EC = (b - a) .* ((X_TestIRI_EC - minX_TrainIRI_EC) ./ (maxX_TrainIRI_EC - minX_TrainIRI_EC)) + a;

% Specify N for CV folds
N_Folds = 10;

% Specify 1,000 lambdas
lambdas = logspace(-3, 0, 1000);

% Fit ridge regression
X_TrainTrans = X_TrainIRI_EC';
rng default
modelIRI_EC = fitrlinear(X_TrainTrans, Y_TrainIRI_EC, 'ObservationsIn', 'columns', 'KFold', N_Folds, 'Lambda', lambdas, 'Learner', 'leastsquares', 'Regularization', 'ridge', 'Solver', 'lbfgs');
MSEsIRI_EC = kfoldLoss(modelIRI_EC);

% Fit the same ridge regression but without CV
modelNonCV_IRI_EC = fitrlinear(X_TrainTrans, Y_TrainIRI_EC, 'ObservationsIn', 'columns', 'Lambda', lambdas, 'Learner', 'leastsquares', 'Regularization', 'ridge', 'Solver', 'lbfgs');

% Extract the non-CV betas
modelNonCV_BetasIRI_EC = modelNonCV_IRI_EC.Beta;

% Select the modelIRI_EC corresponding to the minimal MSE
[MSE_MinIRI_EC, MSE_MinID_IRI_EC] = min(MSEsIRI_EC);
modelFinalIRI_EC = selectModels(modelNonCV_IRI_EC, MSE_MinID_IRI_EC);

% Extract the final betas and intercept
modelFinalBetasIRI_EC = modelFinalIRI_EC.Beta;
modelFinalInterceptIRI_EC = modelFinalIRI_EC.Bias;

% Predict Y in the test set
Y_TestHatIRI_EC = cat(2, ones(size(X_TestIRI_EC, 1), 1), X_TestIRI_EC) * [modelFinalInterceptIRI_EC; modelFinalBetasIRI_EC];

% Calculate stats
L2_NormIRI_EC = norm(modelFinalBetasIRI_EC); % L2 norm
SSE_IRI_EC = sum((Y_TestIRI_EC - Y_TestHatIRI_EC) .^ 2); % Sum of Squared Errors
SST_IRI_EC = sum((Y_TestIRI_EC - mean(Y_TestIRI_EC)) .^ 2); % Total Sum of Squares
MSE_IRI_EC = mean((Y_TestIRI_EC - Y_TestHatIRI_EC) .^ 2); % Mean Squared Error
MAE_IRI_EC = mean(abs(Y_TestIRI_EC - Y_TestHatIRI_EC)); % Mean Absolute Error
R_SquaredIRI_EC = 1 - (SSE_IRI_EC / SST_IRI_EC) % Coefficient of Determination

% save ridgeIRI_EC_SA

%% 2.3 Fit the model in the training set: PCL-R F1

clc; clear; close all

load dataTrainTestSA

% Specify X and Y data
X_TrainPCL_Factor1 = X_Train;
X_TestPCL_Factor1 = X_Test;
Y_PCL_Factor1 = data.PCL_Factor1;

% Partition Y data
Y_TrainPCL_Factor1 = Y_PCL_Factor1(CV.training, :);
Y_TestPCL_Factor1 = Y_PCL_Factor1(CV.test, :);

clearvars -except *PCL_Factor1

% Normalise X data
a = 0;
b = 1;

minX_TrainPCL_Factor1 = min(X_TrainPCL_Factor1);
maxX_TrainPCL_Factor1 = max(X_TrainPCL_Factor1);

X_TrainPCL_Factor1 = (b - a) .* ((X_TrainPCL_Factor1 - minX_TrainPCL_Factor1) ./ (maxX_TrainPCL_Factor1 - minX_TrainPCL_Factor1)) + a;
X_TestPCL_Factor1 = (b - a) .* ((X_TestPCL_Factor1 - minX_TrainPCL_Factor1) ./ (maxX_TrainPCL_Factor1 - minX_TrainPCL_Factor1)) + a;

% Specify N for CV folds
N_Folds = 10;

% Specify 1,000 lambdas
lambdas = logspace(-3, 0, 1000);

% Fit ridge regression
X_TrainTrans = X_TrainPCL_Factor1';
rng default
modelPCL_Factor1 = fitrlinear(X_TrainTrans, Y_TrainPCL_Factor1, 'ObservationsIn', 'columns', 'KFold', N_Folds, 'Lambda', lambdas, 'Learner', 'leastsquares', 'Regularization', 'ridge', 'Solver', 'lbfgs');
MSEsPCL_Factor1 = kfoldLoss(modelPCL_Factor1);

% Fit the same ridge regression but without CV
modelNonCV_PCL_Factor1 = fitrlinear(X_TrainTrans, Y_TrainPCL_Factor1, 'ObservationsIn', 'columns', 'Lambda', lambdas, 'Learner', 'leastsquares', 'Regularization', 'ridge', 'Solver', 'lbfgs');

% Extract the non-CV betas
modelNonCV_BetasPCL_Factor1 = modelNonCV_PCL_Factor1.Beta;

% Select the modelPCL_Factor1 corresponding to the minimal MSE
[MSE_MinPCL_Factor1, MSE_MinID_PCL_Factor1] = min(MSEsPCL_Factor1);
modelFinalPCL_Factor1 = selectModels(modelNonCV_PCL_Factor1, MSE_MinID_PCL_Factor1);

% Extract the final betas and intercept
modelFinalBetasPCL_Factor1 = modelFinalPCL_Factor1.Beta;
modelFinalInterceptPCL_Factor1 = modelFinalPCL_Factor1.Bias;

% Predict Y in the test set
Y_TestHatPCL_Factor1 = cat(2, ones(size(X_TestPCL_Factor1, 1), 1), X_TestPCL_Factor1) * [modelFinalInterceptPCL_Factor1; modelFinalBetasPCL_Factor1];

% Calculate stats
L2_NormPCL_Factor1 = norm(modelFinalBetasPCL_Factor1); % L2 norm
SSE_PCL_Factor1 = sum((Y_TestPCL_Factor1 - Y_TestHatPCL_Factor1) .^ 2); % Sum of Squared Errors
SST_PCL_Factor1 = sum((Y_TestPCL_Factor1 - mean(Y_TestPCL_Factor1)) .^ 2); % Total Sum of Squares
MSE_PCL_Factor1 = mean((Y_TestPCL_Factor1 - Y_TestHatPCL_Factor1) .^ 2); % Mean Squared Error
MAE_PCL_Factor1 = mean(abs(Y_TestPCL_Factor1 - Y_TestHatPCL_Factor1)); % Mean Absolute Error
R_SquaredPCL_Factor1 = 1 - (SSE_PCL_Factor1 / SST_PCL_Factor1); % Coefficient of Determination

% save ridgePCL_Factor1_SA

%% 2.4 Fit the model in the training set: PCL-R F2

clc; clear; close all

load dataTrainTestPCL_Factor2_SA

% Specify X and Y data
X_TrainPCL_Factor2 = X_Train;
X_TestPCL_Factor2 = X_Test;
Y_PCL_Factor2 = data.PCL_Factor2;

% Partition Y data
Y_TrainPCL_Factor2 = Y_PCL_Factor2(CV.training, :);
Y_TestPCL_Factor2 = Y_PCL_Factor2(CV.test, :);

clearvars -except *PCL_Factor2

% Normalise X data
a = 0;
b = 1;

minX_TrainPCL_Factor2 = min(X_TrainPCL_Factor2);
maxX_TrainPCL_Factor2 = max(X_TrainPCL_Factor2);

X_TrainPCL_Factor2 = (b - a) .* ((X_TrainPCL_Factor2 - minX_TrainPCL_Factor2) ./ (maxX_TrainPCL_Factor2 - minX_TrainPCL_Factor2)) + a;
X_TestPCL_Factor2 = (b - a) .* ((X_TestPCL_Factor2 - minX_TrainPCL_Factor2) ./ (maxX_TrainPCL_Factor2 - minX_TrainPCL_Factor2)) + a;

% Specify N for CV folds
N_Folds = 10;

% Specify 1,000 lambdas
lambdas = logspace(-3, 0, 1000);

% Fit ridge regression
X_TrainTrans = X_TrainPCL_Factor2';
rng default
modelPCL_Factor2 = fitrlinear(X_TrainTrans, Y_TrainPCL_Factor2, 'ObservationsIn', 'columns', 'KFold', N_Folds, 'Lambda', lambdas, 'Learner', 'leastsquares', 'Regularization', 'ridge', 'Solver', 'lbfgs');
MSEsPCL_Factor2 = kfoldLoss(modelPCL_Factor2);

% Fit the same ridge regression but without CV
modelNonCV_PCL_Factor2 = fitrlinear(X_TrainTrans, Y_TrainPCL_Factor2, 'ObservationsIn', 'columns', 'Lambda', lambdas, 'Learner', 'leastsquares', 'Regularization', 'ridge', 'Solver', 'lbfgs');

% Extract the non-CV betas
modelNonCV_BetasPCL_Factor2 = modelNonCV_PCL_Factor2.Beta;

% Select the modelPCL_Factor2 corresponding to the minimal MSE
[MSE_MinPCL_Factor2, MSE_MinID_PCL_Factor2] = min(MSEsPCL_Factor2);
modelFinalPCL_Factor2 = selectModels(modelNonCV_PCL_Factor2, MSE_MinID_PCL_Factor2);

% Extract the final betas and intercept
modelFinalBetasPCL_Factor2 = modelFinalPCL_Factor2.Beta;
modelFinalInterceptPCL_Factor2 = modelFinalPCL_Factor2.Bias;

% Predict Y in the test set
Y_TestHatPCL_Factor2 = cat(2, ones(size(X_TestPCL_Factor2, 1), 1), X_TestPCL_Factor2) * [modelFinalInterceptPCL_Factor2; modelFinalBetasPCL_Factor2];

% Calculate stats
L2_NormPCL_Factor2 = norm(modelFinalBetasPCL_Factor2); % L2 norm
SSE_PCL_Factor2 = sum((Y_TestPCL_Factor2 - Y_TestHatPCL_Factor2) .^ 2); % Sum of Squared Errors
SST_PCL_Factor2 = sum((Y_TestPCL_Factor2 - mean(Y_TestPCL_Factor2)) .^ 2); % Total Sum of Squares
MSE_PCL_Factor2 = mean((Y_TestPCL_Factor2 - Y_TestHatPCL_Factor2) .^ 2); % Mean Squared Error
MAE_PCL_Factor2 = mean(abs(Y_TestPCL_Factor2 - Y_TestHatPCL_Factor2)); % Mean Absolute Error
R_SquaredPCL_Factor2 = 1 - (SSE_PCL_Factor2 / SST_PCL_Factor2) % Coefficient of Determination

% save ridgePCL_Factor2_SA

%% 3.1 Evaluate the model in the test set: IRI-PT

clc; clear; close all

load ridgeIRI_PT_SA

N_Perms = 10000;
Y_TestHatPermIRI_PT = nan(size(Y_TestHatIRI_PT, 1), N_Perms);
MSE_PermIRI_PT = nan(N_Perms, 1);

rng default
for i = 1 : N_Perms
    
    disp(i)
    X_TestPermIRI_PT = X_TestIRI_PT(randperm(size(X_TestIRI_PT, 1)), :);
    Y_TestHatPermIRI_PT(:, i) = cat(2, ones(size(X_TestIRI_PT, 1), 1), X_TestPermIRI_PT) * [modelFinalInterceptIRI_PT; modelFinalBetasIRI_PT];
    MSE_PermIRI_PT(i) = mean((Y_TestIRI_PT - Y_TestHatPermIRI_PT(:, i)) .^ 2);
    
end

MSE_DistroIRI_PT = [MSE_IRI_PT; MSE_PermIRI_PT];
MSE_DistroRankIRI_PT = tiedrank(MSE_DistroIRI_PT) / size(MSE_DistroIRI_PT, 1); % One-tailed
P_PermIRI_PT = MSE_DistroRankIRI_PT(1)

% save ridgeIRI_PT_SA_Perm

%% 3.2 Evaluate the model in the test set: PCL-R F1

clc; clear; close all

load ridgePCL_Factor1_SA

N_Perms = 10000;
Y_TestHatPermPCL_Factor1 = nan(size(Y_TestHatPCL_Factor1, 1), N_Perms);
MSE_PermPCL_Factor1 = nan(N_Perms, 1);

rng default
for i = 1 : N_Perms
    
    disp(i)
    X_TestPermPCL_Factor1 = X_TestPCL_Factor1(randperm(size(X_TestPCL_Factor1, 1)), :);
    Y_TestHatPermPCL_Factor1(:, i) = cat(2, ones(size(X_TestPCL_Factor1, 1), 1), X_TestPermPCL_Factor1) * [modelFinalInterceptPCL_Factor1; modelFinalBetasPCL_Factor1];
    MSE_PermPCL_Factor1(i) = mean((Y_TestPCL_Factor1 - Y_TestHatPermPCL_Factor1(:, i)) .^ 2);
    
end

MSE_DistroPCL_Factor1 = [MSE_PCL_Factor1; MSE_PermPCL_Factor1];
MSE_DistroRankPCL_Factor1 = tiedrank(MSE_DistroPCL_Factor1) / size(MSE_DistroPCL_Factor1, 1); % One-tailed
P_PermPCL_Factor1 = MSE_DistroRankPCL_Factor1(1)

% save ridgePCL_Factor1_SA_Perm

% Calculate a 95% CI for R2
N_Boot = 10000; % N for bootstrap samples
bootR2 = zeros(N_Boot, 1); % To store R^2 for each bootstrap sample
N_Test = height(Y_TestPCL_Factor1); % N for test samples

rng default
for i = 1 : N_Boot
    
    % Generate bootstrap sample indices (with replacement)
    bootIDX = randsample(N_Test, N_Test, true);
    
    % Extract bootstrap test data
    X_TestBoot = X_TestPCL_Factor1(bootIDX, :);
    Y_TestBoot = Y_TestPCL_Factor1(bootIDX);
    
    % Predict Y for bootstrap test set
    Y_TestHatBoot = cat(2, ones(size(X_TestBoot, 1), 1), X_TestBoot) * ...
        [modelFinalInterceptPCL_Factor1; modelFinalBetasPCL_Factor1];
    
    % Calculate SSE and SST for bootstrap sample
    SSE_Boot = sum((Y_TestBoot - Y_TestHatBoot) .^ 2);
    SST_Boot = sum((Y_TestBoot - mean(Y_TestBoot)) .^ 2);
    
    % Calculate R^2 for bootstrap sample
    bootR2(i) = 1 - (SSE_Boot / SST_Boot);

end

% Compute 95% confidence interval using percentile method
CI_R2 = prctile(bootR2, [2.5, 97.5]);

% Display results
fprintf('Out-of-sample R^2: %.2f\n', R_SquaredPCL_Factor1);
fprintf('95%% Confidence Interval for R^2: [%.2f, %.2f]\n', CI_R2(1), CI_R2(2));