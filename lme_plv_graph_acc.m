function lme_plv_graph_acc(freq)
% LME_PLV_GRAPH_ACC - linear mixed-effects model with dependent variable: 
% difference in PLV node strength between stimulus types; fixed effects: 
% correct/error trial, Yeo-7 network (electrode localization); random 
% effects: subject, electrode nested in subject.
%
% Yeo-7 network localization data must be accessed from OSF:
% https://doi.org/10.17605/OSF.IO/RX2ZD
%
% Inputs:
% freq = 'theta'/'alpha'/'beta' PLV node strength frequency
%
% Example:
% lme_plv_graph_acc('theta')
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

    % load HFB data for correct and error trials, in try/catch in case of
    % no error trials
    if isfile(fullfile(datdir, 'plv', ['cf_' freq '_graph.mat']))
        cf_cor = load(fullfile(datdir, 'plv', ['cf_' freq '_graph'], freq));
        cf_inc = load(fullfile(datdir, 'plv', ['cf_' freq '_graph_inc'], freq));

        cl_cor = load(fullfile(datdir, 'plv', ['cl_' freq '_graph'], freq));
        try
            cl_inc = load(fullfile(datdir, 'plv', ['cl_' freq '_graph_inc'], freq));
        catch
        end
        
        ch_idx = ismember(ch, cf_cor.(freq).label);
        yeo_key = yeo_key(ch_idx);      
        
        for e = 1:length(cf_cor.(freq).label)
            c = c + 1;
            sub{c} = sbj{s};
            elec{c} = cf_cor.(freq).label{e};
            network{c} = yeo_key{e};
            cf_cor_m(c) = cf_cor.(freq).diff_sum(e);
            cf_inc_m(c) = cf_inc.(freq).diff_sum(e);
            cl_cor_m(c) = cl_cor.(freq).diff_sum(e);
            if ~strcmp(sbj{s},'S3')
                cl_inc_m(c) = cl_inc.(freq).diff_sum(e);
            else
                cl_inc_m(c) = nan;
            end
        end

        clear yeo_key *cor *inc datdir ch_idx
    end
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
save(fullfile(savdir, ['lme_' freq '_graph_acc']), 'model_data', '*lme*');

end
