%% Figure 6

%% 1.1 Test for a median class difference: CT

clc; clear; close all

load('gradientsCT_Low', 'G1_CT_Low')
load('gradientsCT_High', 'G1_CT_High')
load classesMesulamGlasser

clearvars -except G1_CT_High G1_CT_Low Mesulam*

MesulamLabels = {'P' 'H' 'U' 'I'};

N_Tests = max(MesulamGlasser);
for i = 1 : max(MesulamGlasser)
    
    G1_CT_MesulamName = sprintf('G1_CT_Mesulam%.0f', i);
    
    [P_Value, ~, stats] = signrank(G1_CT_High(MesulamGlasser == i), G1_CT_Low(MesulamGlasser == i));
    
    diff = G1_CT_High(MesulamGlasser == i) - G1_CT_Low(MesulamGlasser == i);
    N = sum(diff ~= 0 & ~isnan(diff));
    G1_CT_EffectSizeR = stats.zval / sqrt(N);
    
    eval([G1_CT_MesulamName ' = G1_CT_EffectSizeR;']);
    
    P_ValueBon = P_Value * N_Tests;
    if P_ValueBon > 1
    P_ValueBon = 1;
    end
    
    G1_CT_P_ValueBonMesulamName = sprintf('G1_CT_P_ValueBonMesulam%.0f', i);
    eval([G1_CT_P_ValueBonMesulamName ' = P_ValueBon;'])
    
    eval(sprintf('G1_CT_P_ValueBonMesulam%.0f', i))

end

clearvars -except G1* Mesulam* MesulamLabels

% save G1_CT_MesulamEffectSizeR

%% 1.2 Test for a median network difference: CT

clc; clear; close all

load('gradientsCT_Low', 'G1_CT_Low')
load('gradientsCT_High', 'G1_CT_High')
load networksYeoGlasser

clearvars -except G1_CT_High G1_CT_Low Yeo*

YeoLabels = {'V' 'S' 'DA' 'VA' 'L' 'F' 'DM'};

N_Tests = max(YeoGlasser);
for i = 1 : max(YeoGlasser)
    
    G1_CT_YeoName = sprintf('G1_CT_Yeo%.0f', i);
    
    [P_Value, ~, stats] = signrank(G1_CT_High(YeoGlasser == i), G1_CT_Low(YeoGlasser == i));
    
    diff = G1_CT_High(YeoGlasser == i) - G1_CT_Low(YeoGlasser == i);
    N = sum(diff ~= 0 & ~isnan(diff));
    G1_CT_EffectSizeR = stats.zval / sqrt(N);
    
    eval([G1_CT_YeoName ' = G1_CT_EffectSizeR;']);
    
    P_ValueBon = P_Value * N_Tests;
    if P_ValueBon > 1
    P_ValueBon = 1;
    end
    
    G1_CT_P_ValueBonYeoName = sprintf('G1_CT_P_ValueBonYeo%.0f', i);
    eval([G1_CT_P_ValueBonYeoName ' = P_ValueBon;'])
    
    eval(sprintf('G1_CT_P_ValueBonYeo%.0f', i))

end

clearvars -except G1* Yeo* YeoLabels

% save G1_CT_YeoEffectSizeR

%% 1.3 Test for a median class difference: SA

clc; clear; close all

load('gradientsSA_Low', 'G1_SA_Low')
load('gradientsSA_High', 'G1_SA_High')
load classesMesulamGlasser

clearvars -except G1_SA_High G1_SA_Low Mesulam*

MesulamLabels = {'P' 'H' 'U' 'I'};

N_Tests = max(MesulamGlasser);
for i = 1 : max(MesulamGlasser)
    
    G1_SA_MesulamName = sprintf('G1_SA_Mesulam%.0f', i);
    
    [P_Value, ~, stats] = signrank(G1_SA_High(MesulamGlasser == i), G1_SA_Low(MesulamGlasser == i));
    
    diff = G1_SA_High(MesulamGlasser == i) - G1_SA_Low(MesulamGlasser == i);
    N = sum(diff ~= 0 & ~isnan(diff));
    G1_SA_EffectSizeR = stats.zval / sqrt(N);
    
    eval([G1_SA_MesulamName ' = G1_SA_EffectSizeR;']);
    
    P_ValueBon = P_Value * N_Tests;
    if P_ValueBon > 1
    P_ValueBon = 1;
    end
    
    G1_SA_P_ValueBonMesulamName = sprintf('G1_SA_P_ValueBonMesulam%.0f', i);
    eval([G1_SA_P_ValueBonMesulamName ' = P_ValueBon;'])
    
    eval(sprintf('G1_SA_P_ValueBonMesulam%.0f', i))

end

clearvars -except G1* Mesulam* MesulamLabels

% save G1_SA_MesulamEffectSizeR

%% 1.4 Test for a median network difference: SA

clc; clear; close all

load('gradientsSA_Low', 'G1_SA_Low')
load('gradientsSA_High', 'G1_SA_High')
load networksYeoGlasser

clearvars -except G1_SA_High G1_SA_Low Yeo*

YeoLabels = {'V' 'S' 'DA' 'VA' 'L' 'F' 'DM'};

N_Tests = max(YeoGlasser);
for i = 1 : max(YeoGlasser)
    
    G1_SA_YeoName = sprintf('G1_SA_Yeo%.0f', i);
    
    [P_Value, ~, stats] = signrank(G1_SA_High(YeoGlasser == i), G1_SA_Low(YeoGlasser == i));
    
    diff = G1_SA_High(YeoGlasser == i) - G1_SA_Low(YeoGlasser == i);
    N = sum(diff ~= 0 & ~isnan(diff));
    G1_SA_EffectSizeR = stats.zval / sqrt(N);
    
    eval([G1_SA_YeoName ' = G1_SA_EffectSizeR;']);
    
    P_ValueBon = P_Value * N_Tests;
    if P_ValueBon > 1
    P_ValueBon = 1;
    end
    
    G1_SA_P_ValueBonYeoName = sprintf('G1_SA_P_ValueBonYeo%.0f', i);
    eval([G1_SA_P_ValueBonYeoName ' = P_ValueBon;'])
    
    eval(sprintf('G1_SA_P_ValueBonYeo%.0f', i))

end

clearvars -except G1* Yeo* YeoLabels

% save G1_SA_YeoEffectSizeR