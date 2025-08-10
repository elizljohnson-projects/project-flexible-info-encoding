% Linear mixed-effects model with dependent variable: difference in HFB
% power between stimulus types; fixed effects: correct/error trial, Yeo-7 
% network (electrode localization); random effects: subject, electrode
% nested in subject.
%
% Yeo-7 network localization data must be accessed from OSF:
% https://doi.org/10.17605/OSF.IO/RX2ZD
%
% Copyright (c) 2023
% EL Johnson, PhD

% set directories
pth = pwd;
savdir = fullfile(pth, 'hfb');
mkdir(savdir);

% subject list
sbj = {'S1','S2','S3','S4','S5','S6','S7','S8','S9','S10','S11'};

% construct model data table with data from all subjects
c = 0;
for s = 1:length(sbj)
    datdir = fullfile(pth, 'data', sbj{s});

    % load Yeo-7 network data
    load(fullfile(datdir, 'yeo'), 'yeo_key');

    % load HFB data for correct and error trials, in try/catch in case of
    % no error trials
    cf_cor = load(fullfile(datdir, 'hfb', 'cf_hfb'), 'hfb');
    cf_inc = load(fullfile(datdir, 'hfb', 'cf_hfb_inc'), 'hfb');

    cl_cor = load(fullfile(datdir, 'hfb', 'cl_hfb'), 'hfb');
    try
        cl_inc = load(fullfile(datdir, 'hfb', 'cl_hfb_inc'), 'hfb');
    catch
    end
    
    for e = 1:length(cf_cor.hfb.label)
        c = c + 1;
        sub{c} = sbj{s};
        elec{c} = cf_cor.hfb.label{e};
        network{c} = yeo_key{e};
        cf_cor_m(c) = cf_cor.hfb.diff(e);
        cf_inc_m(c) = cf_inc.hfb.diff(e);
        cl_cor_m(c) = cl_cor.hfb.diff(e);
        try
            cl_inc_m(c) = cl_inc.hfb.diff(e);
        catch
            cl_inc_m(c) = nan;
        end
    end
    
    clear yeo_key *cor *inc datdir
end

acc = [zeros(length(sub),1); ones(length(sub),1)];
cf_m = [cf_inc_m'; cf_cor_m'];
cl_m = [cl_inc_m'; cl_cor_m'];
sub = [sub'; sub'];
elec = [elec'; elec'];
network = [network'; network'];

model_data = table(sub, elec, network, cf_m, cl_m, acc, ...
    'VariableNames', {'sub', 'elec', 'network', 'cf_m', 'cl_m', 'acc'});

% run model on context first condition
cf_lme = fitlme(model_data, 'cf_m ~ acc * network + (1|sub) + (1|sub:elec)');
cf_lme_anova = anova(cf_lme);

% run model on context last condition
cl_lme = fitlme(model_data, 'cl_m ~ acc * network + (1|sub) + (1|sub:elec)');
cl_lme_anova = anova(cl_lme);

% save model data table and results
save(fullfile(savdir, 'lme_acc'), 'model_data', '*lme*');
