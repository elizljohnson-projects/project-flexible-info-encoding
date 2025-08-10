function lme_plv_graph_acc_int(freq)
% LME_PLV_GRAPH_ACC_INT - linear mixed-effects model for post hoc testing 
% the accuracy * network interaction on the difference in PLV node strength
% between stimulus types (i.e., results of lme_plv_graph_acc). Run the 
% model on each network with fixed effects: correct/error trial; random 
% effects: subject, electrode nested in subject.
%
% Inputs:
% freq = 'theta'/'alpha'/'beta' PLV node strength frequency
%
% Example:
% lme_plv_graph_acc_int('theta')
%
% Copyright (c) 2023
% EL Johnson, PhD

clearvars -except freq

% set directories
pth = pwd;
datdir = fullfile(pth, 'plv'); % output of lme_plv_graph_acc
savdir = datdir;

% load model data table
load(fullfile(datdir, ['lme_' freq '_graph_acc']), 'model_data');

% index each network in the table
network{1} = strcmp(model_data.network, 'VIS');
network{2} = strcmp(model_data.network, 'SM');
network{3} = strcmp(model_data.network, 'DAN');
network{4} = strcmp(model_data.network, 'VAN');
network{5} = strcmp(model_data.network, 'LB');
network{6} = strcmp(model_data.network, 'FP');
network{7} = strcmp(model_data.network, 'DM');

% initialize variables for model outputs
cf_lme_acc = cell(length(network), 1);
cf_lme_acc_anova = cf_lme_acc;
cf_lme_acc_t = nan(length(network), 1);

cl_lme_acc = cf_lme_acc;
cl_lme_acc_anova = cf_lme_acc;
cl_lme_acc_t = cf_lme_acc_t;

lme_acc_network = cf_lme_acc;

% loop through networks
for n1 = 1:length(network)
    tmp_model = model_data(network{n1},:);
    
    lme_acc_network{n1} = num2str(n1);
    
    % run model on context first condition, in try/catch in case of
    % missing data from no electrodes in network
    try
        cf_lme_acc{n1} = fitlme(tmp_model, 'cf_m ~ acc + (1|sub) + (1|sub:elec)');
        cf_lme_acc_anova{n1} = anova(cf_lme_acc{n1});
        cf_lme_acc_t(n1) = cf_lme_acc{n1}.Coefficients.tStat(2);
    catch
        cf_lme_acc{n1} = nan;
        cf_lme_acc_anova{n1} = nan;
        cf_lme_acc_t(n1) = nan;
    end
    
    % run model on context first condition, in try/catch in case of
    % missing data from no electrodes in network
    try
        cl_lme_acc{n1} = fitlme(tmp_model, 'cl_m ~ acc + (1|sub) + (1|sub:elec)');
        cl_lme_acc_anova{n1} = anova(cl_lme_acc{n1});
        cl_lme_acc_t(n1) = cl_lme_acc{n1}.Coefficients.tStat(2);
    catch
        cl_lme_acc{n1} = nan;
        cl_lme_acc_anova{n1} = nan;
        cl_lme_acc_t(n1) = nan;
    end
    
    clear tmp*
end

% save model results
save(fullfile(savdir, ['lme_' freq '_graph_acc_int']), '*lme*');

end
