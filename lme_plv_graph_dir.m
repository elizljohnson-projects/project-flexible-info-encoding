function lme_plv_graph_dir(freq)
% LME_PLV_GRAPH_DIR - linear mixed-effects model with dependent variable: 
% correlation direction of PLV node strength difference with RT; fixed 
% effects: correlation across trials, Yeo-7 network (electrode 
% localization); random effects: subject, electrode nested in subject.
%
% Yeo-7 network localization data must be accessed from OSF:
% https://doi.org/10.17605/OSF.IO/RX2ZD
%
% Inputs:
% freq = 'theta'/'alpha'/'beta' PLV node strength frequency
%
% Example:
% lme_plv_graph_dir('theta')
%
% Copyright (c) 2023
% EL Johnson, PhD

clearvars -except freq

% set directories
pth = pwd;
savdir = fullfile(pth, 'plv');
mkdir(savdir);

% subject list
sbj = {'S1','S2','S3','S4','S5','S6','S7','S8','S9','S10','S11'};

% construct model data table with data from all subjects
c = 0;
for s = 1:length(sbj)
    datdir = fullfile(pth, 'data', sbj{s});

    % load Yeo-7 network data
    load(fullfile(datdir, 'yeo'), 'yeo_key');

    % load channel key for Yeo-7 list
    hfb = load(fullfile(datdir, 'hfb', 'cf_hfb'));
    ch = hfb.hfb.label; clear hfb
    
    if isfile(fullfile(datdir, 'plv', ['cf_' freq '_graph.mat']))
        cf = load(fullfile(datdir, 'plv', ['cf_' freq '_graph']), freq);
        cl = load(fullfile(datdir, 'cl_', [freq '_graph']), freq);
        
        ch_idx = ismember(ch, cf.(freq).label);
        yeo_key = yeo_key(ch_idx);

        for e = 1:length(cf.(freq).label)
            c = c + 1;
            sub{c} = sbj{s};
            elec{c} = cf.(freq).label{e};
            network{c} = yeo_key{e};
            cf_z(c) = cf.(freq).r_sum(e);
            cl_z(c) = cl.(freq).r_sum(e);
        end
    end
    
    clear yeo_key cf cl datdir ch_idx
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
save(fullfile(savdir, ['lme_' freq '_graph_dir']), 'model_data', '*lme*');

end
