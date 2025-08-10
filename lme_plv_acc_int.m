function lme_plv_acc_int(freq)
% LME_PLV_ACC_INT - linear mixed-effects model for post hoc testing the
% accuracy * network interaction on the difference in PLV between stimulus
% types (i.e., results of lme_plv_acc). Run the model on each network with
% fixed effects: correct/error trial; random effects: subject, electrode
% nested in subject.
%
% Inputs:
% freq = 'theta'/'alpha'/'beta' PLV frequency
%
% Example:
% lme_plv_acc_int('theta')
%
% Copyright (c) 2023
% EL Johnson, PhD

clearvars -except freq

% set directories
pth = pwd;
datdir = fullfile(pth, 'plv'); % output of lme_plv_acc
savdir = datdir;

% load model data table
load(fullfile(datdir, ['lme_' freq '_acc']), 'model_data');

% index each network in the table
network1{1} = strcmp(model_data.network1, 'VIS');
network1{2} = strcmp(model_data.network1, 'SM');
network1{3} = strcmp(model_data.network1, 'DAN');
network1{4} = strcmp(model_data.network1, 'VAN');
network1{5} = strcmp(model_data.network1, 'LB');
network1{6} = strcmp(model_data.network1, 'FP');
network1{7} = strcmp(model_data.network1, 'DM');

network2{1} = strcmp(model_data.network2, 'VIS');
network2{2} = strcmp(model_data.network2, 'SM');
network2{3} = strcmp(model_data.network2, 'DAN');
network2{4} = strcmp(model_data.network2, 'VAN');
network2{5} = strcmp(model_data.network2, 'LB');
network2{6} = strcmp(model_data.network2, 'FP');
network2{7} = strcmp(model_data.network2, 'DM');

% initialize variables for model outputs
cf_lme_acc = cell(length(network1), length(network2));
cf_lme_acc_anova = cf_lme_acc;
cf_lme_acc_t = nan(length(network1), length(network2));

cl_lme_acc = cf_lme_acc;
cl_lme_acc_anova = cf_lme_acc;
cl_lme_acc_t = cf_lme_acc_t;

lme_acc_network = cf_lme_acc;

% loop through networks
for n1 = 1:length(network1)
    for n2 = 1:length(network2)
        if n1 >= n2
            tmp_idx = network1{n1} + network2{n2} == 2;
            tmp_model = model_data(tmp_idx,:);
            
            lme_acc_network{n1,n2} = strcat(num2str(n1), '-', num2str(n2));
            
            % run model on context first condition, in try/catch in case of
            % missing data from no electrodes in network
            try
                cf_lme_acc{n1,n2} = fitlme(tmp_model, ...
                    'cf_m ~ acc + (1|sub) + (1|sub:elec1) + (1|sub:elec2)');
                cf_lme_acc_anova{n1,n2} = anova(cf_lme_acc{n1,n2});
                cf_lme_acc_t(n1,n2) = cf_lme_acc{n1,n2}.Coefficients.tStat(2);
            catch
                cf_lme_acc{n1,n2} = nan;
                cf_lme_acc_anova{n1,n2} = nan;
                cf_lme_acc_t(n1,n2) = nan;
            end
            
            % run model on context first condition, in try/catch in case of
            % missing data from no electrodes in network
            try
                cl_lme_acc{n1,n2} = fitlme(tmp_model, ...
                    'cl_m ~ acc + (1|sub) + (1|sub:elec1) + (1|sub:elec2)');
                cl_lme_acc_anova{n1,n2} = anova(cl_lme_acc{n1,n2});
                cl_lme_acc_t(n1,n2) = cl_lme_acc{n1,n2}.Coefficients.tStat(2);
            catch
                cl_lme_acc{n1,n2} = nan;
                cl_lme_acc_anova{n1,n2} = nan;
                cl_lme_acc_t(n1,n2) = nan;
            end

            clear tmp*
        end
    end
end

% save model results
save(fullfile(savdir, ['lme_' freq '_acc_int']), '*lme*');

end
