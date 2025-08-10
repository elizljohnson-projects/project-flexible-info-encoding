% Linear mixed-effects model with dependent variable: correlation direction
% of HFB power difference with RT; fixed effects: correlation across 
% trials, Yeo-7 network (electrode localization); random effects: subject, 
% electrode nested in subject.
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

    % load HFB correlation across trials
    cf = load(fullfile(datdir, 'hfb', 'cf_hfb'), 'hfb');
    cl = load(fullfile(datdir, 'hfb', 'cl_hfb'), 'hfb');
    
    for e = 1:length(cf.hfb.label)
        c = c + 1;
        sub{c} = sbj{s};
        elec{c} = cf.hfb.label{e};
        network{c} = yeo_key{e};
        cf_z(c) = cf.hfb.r_z(e);
        cl_z(c) = cl.hfb.r_z(e);
    end
    
    clear yeo_key cf cl datdir
end

dir = [zeros(length(sub),1); ones(length(sub),1)];
cf_z = [zeros(length(sub),1); cf_z'];
cl_z = [zeros(length(sub),1); cl_z'];
sub = [sub'; sub'];
elec = [elec'; elec'];
network = [network'; network'];

model_data = table(sub, elec, network, cf_z, cl_z, dir, ...
    'VariableNames', {'sub', 'elec', 'network', 'cf_z', 'cl_z', 'dir'});

% run model on context first condition
cf_lme = fitlme(model_data, 'cf_z ~ dir * network + (1|sub) + (1|sub:elec)');
cf_lme_anova = anova(cf_lme);

% run model on context last condition
cl_lme = fitlme(model_data, 'cl_z ~ dir * network + (1|sub) + (1|sub:elec)');
cl_lme_anova = anova(cl_lme);

% save model data table and results
save(fullfile(savdir, 'lme_dir'), 'model_data', '*lme*');
