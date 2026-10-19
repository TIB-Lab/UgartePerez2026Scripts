clear all
close all

%% Load data
T = readtable('data.xls','VariableNamingRule','preserve');    

y = str2double(T{:,1});

animal = string(T{:,2});     

region_full = string(T{:,3});        

%% Extract SEX
sex = strings(size(animal));
second_letter = extractBetween(animal,2,2); % get 2nd character

sex(second_letter == "F") = "F";
sex(second_letter == "M") = "M";
sex = categorical(sex);


%% Extract BATCH
batch = strings(size(animal));
first_letter = extractBetween(animal,1,1); % get 1st character

batch(first_letter == "E") = "2";
batch(first_letter ~= "E") = "1";
batch = categorical(batch);


%% Extract ANIMAL ID
animal = categorical(animal);
animal_clean=animal(keep);

%% Extract REGION (Prelim a / Prelim b / RSA)
region = strings(size(region_full));

region(contains(region_full,"Prelim") & endsWith(region_full,"a")) = "Prelim_a"; %a
region(contains(region_full,"Prelim") & endsWith(region_full,"b")) = "Prelim_b"; %b
region(contains(region_full,"RSA")) = "RSA";

region = categorical(region);
roi = categorical(region_full);

%% Build table
data = table(y, animal, sex, region, batch, roi, ...
    'VariableNames', {'Value','Animal','Sex','Region','Batch','ROI'});

%% Aggregate per animal × region
G = groupsummary(data, {'Animal','Sex','Region','Batch'}, 'median', 'Value');

lme = fitlme(G, 'median_Value ~ Sex * Region + Batch + (1|Animal)');

% ANOVA
a = anova(lme)

regions = categories(G.Region);
sexes   = categories(G.Sex);

%% Plot

figure; hold on

offset = [-0.15 0.15]; % for F vs M separation

for s = 1:numel(sexes)
    for r = 1:numel(regions)
        
        idx = G.Sex == sexes{s} & G.Region == regions{r};
        yvals = G.median_Value(idx);
        
        x = r + offset(s);

        % scatter points
        scatter(repmat(x, size(yvals)), yvals, 50, 'filled')

        hold on

        % ADD LABELS
        animal_names = string(G.Animal(idx));

        for k = 1:numel(yvals)
            text(x + 0.02, yvals(k), animal_names(k), ...
                'FontSize',8, 'Color',[0.3 0.3 0.3])
        end

        % mean ± SEM
        m = mean(yvals, 'omitnan');
        sem = std(yvals,'omitnan')/sqrt(sum(~isnan(yvals)));
        
        errorbar(x, m, sem, 'k', 'LineWidth',1.5)
    end
end

xticks(1:numel(regions))
xticklabels(regions)
xlabel('Region')
ylabel('Mean count per animal')
legend(sexes, 'Location','best')
title('Sex × Region')

box on

%% Get regions
G = groupsummary(data, {'Animal','Sex','Region','Batch'}, 'mean', 'Value');

regions = categories(G.Region);

pvals = nan(numel(regions),1);
effectSize = nan(numel(regions),1);
CI_low = nan(numel(regions),1);
CI_high = nan(numel(regions),1);
meanDiff = nan(numel(regions),1);
meanDiff_CI_low = nan(numel(regions),1);
meanDiff_CI_high = nan(numel(regions),1);

for r = 1:numel(regions)
    
    idx = G.Region == regions{r};
    subdata = G(idx,:);
    
    % LME
    lme_r = fitlme(subdata, 'mean_Value ~ Sex + (1|Animal)');
    
    % ANOVA
    a = anova(lme_r);
    pvals(r) = a.pValue(2);
    
    %% Estimated sex difference and 95% CI
    
    % Fixed-effect coefficients
    coef = lme_r.Coefficients;
    
    % Identify Sex coefficient
    sex_idx = contains(coef.Name, 'Sex');
    
    beta = coef.Estimate(sex_idx);
    beta_CI = coefCI(lme_r);
    
    meanDiff(r) = beta;
    meanDiff_CI_low(r) = beta_CI(sex_idx,1);
    meanDiff_CI_high(r) = beta_CI(sex_idx,2);
    
    %% Cohen's d
    
    % Residual SD
    residualSD = sqrt(lme_r.MSE);
    
    effectSize(r) = beta / residualSD;
    
    % Approximate 95% CI for d
    se_beta = coef.SE(sex_idx);
    
    se_d = se_beta / residualSD;
    
    CI_low(r) = effectSize(r) - 1.96*se_d;
    CI_high(r) = effectSize(r) + 1.96*se_d;
end

%% FDR correction
pvals_fdr = mafdr(pvals, 'BHFDR', true);

%% Results table

results = table( ...
    regions, ...
    meanDiff, ...
    meanDiff_CI_low, ...
    meanDiff_CI_high, ...
    effectSize, ...
    CI_low, ...
    CI_high, ...
    pvals, ...
    pvals_fdr, ...
    'VariableNames', { ...
    'Region', ...
    'MeanDifference', ...
    'MeanDifference_CI_low', ...
    'MeanDifference_CI_high', ...
    'Cohens_d', ...
    'd_CI_low', ...
    'd_CI_high', ...
    'p_uncorrected', ...
    'p_FDR'});

disp(results)