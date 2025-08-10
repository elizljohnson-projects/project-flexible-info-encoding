function peaks_irasa(sid)
% PEAKS_IRASA - compute oscillatory peaks based on stimulus epoch of all 
% correct trials using IRASA method. Detect peaks using statistical
% bootstrapping.
%
% Analysis described in: Wen, H, Liu, Z. Separating fractal and oscillatory 
% components in the power spectrum of neurophysiological signal. Brain 
% Topography 29 (2016). https://doi.org/10.1007/s10548-015-0448-0 
% 
% Ensure FieldTrip is correcty added to the MATLAB path:
%   addpath <path to fieldtrip home directory>
%   ft_defaults
%
% Inputs:
% sid = subject ID (e.g., 'S1')
%
% Example:
% peaks_irasa('S1')
%
% Copyright (c) 2023
% EL Johnson, PhD

clearvars -except sid

% set directories
pth = pwd;
datdir = fullfile(pth, sid);
savdir = fullfile(datdir, 'irasa');
mkdir(savdir);

% set threshold for peak
zthr = 1.96;

% load data
data = load(fullfile(datdir, sid), 'data');
data = data.data;

srate = data.fsample;

% select correct trials
cor = cell2mat(data.trialinfo(:,6));
trialinfo = data.trialinfo;

% select task-responsive elecs
chans_stim = load(fullfile(savdir, 'hfb', 'chans_stim'));
chans_stim = chans_stim.chans_stim;
chans_resp = load(fullfile(savdir, 'hfb', 'chans_resp'));
chans_resp = chans_resp.chans_resp;

cfg = [];
cfg.trials = find(cor == 1);
cfg.channel = unique([chans_stim chans_resp]);

data = ft_selectdata(cfg, data);
data.trialinfo = trialinfo(cfg.trials,:);
clear cor trialinfo chans*

% select stim onset --> resp onset, no buffers
resp = str2double(data.trialinfo(:,18));
resp = resp./srate;

for r = 1:length(data.trial)
    t1 = nearest(data.time{r}, 0);
    t2 = nearest(data.time{r}, resp(r));
    data.trial{r} = data.trial{r}(:,t1:t2);
    data.time{r} = data.time{r}(1,t1:t2);
    clear t1 t2
end

data.trialinfo = 1:size(data.trialinfo,1);
data.trialinfo = data.trialinfo';

% partition the data into ten overlapping sub-segments
[~,idx] = min(cellfun('size', data.trial, 2)); % find shortest trial
w = data.time{idx}(end)-data.time{idx}(1); % window length
cfg = [];
cfg.length = w*.9;
cfg.overlap = 1-((w-cfg.length)/(10-1));

tmp = ft_redefinetrial(cfg, data);

% perform IRASA and regular spectral analysis
cfg = [];
cfg.foilim = [3 30];
cfg.taper = 'hanning';
cfg.pad = 'nextpow2';
cfg.keeptrials = 'yes';

cfg.method = 'irasa';
frac = ft_freqanalysis(cfg, tmp);

cfg.method = 'mtmfft';
orig = ft_freqanalysis(cfg, tmp);
clear tmp*

% average across the sub-segments
tmpf = {};
tmpo = {};
for rpt = unique(frac.trialinfo)'
    cfg = [];
    cfg.trials = find(frac.trialinfo == rpt);
    cfg.avgoverrpt = 'yes';
    
    tmpf{end+1} = ft_selectdata(cfg, frac);
    tmpo{end+1} = ft_selectdata(cfg, orig);
end

frac = ft_appendfreq([], tmpf{:});
orig = ft_appendfreq([], tmpo{:});
clear tmp*

% average across trials
cfg = [];
cfg.avgoverrpt = 'yes';

frac = ft_selectdata(cfg, frac);
orig = ft_selectdata(cfg, orig);

% subtract the fractal component from the power spectrum
cfg = [];
cfg.parameter = 'powspctrm';
cfg.operation = 'x2-x1';
osci = ft_math(cfg, frac, orig);

% detect peaks per electrode
osci.pk.pk = cell(length(osci.label),1);
osci.pk.fr = osci.pk.pk;
osci.pk.wd = osci.pk.fr;
osci.theta = nan(length(osci.label),1);
osci.alpha = osci.theta;
osci.beta = osci.theta;

% z-score using statistical bootstrapping
osci.powspctrm = zcont_5(osci.powspctrm); 

for e = 1:length(osci.label)
    [tmp_1, tmp_fr, tmp_wd, tmp_pk] = findpeaks(osci.powspctrm(e,:), osci.freq);
    tmp1 = tmp_1 > 0; tmp2 = tmp_pk >= zthr; tmp = tmp1 + tmp2 == 2;
    tmp_pk = tmp_pk(tmp); tmp_fr = tmp_fr(tmp); tmp_wd = tmp_wd(tmp);
    [~,tmp] = sort(tmp_pk, 'descend');
    tmp_pk = tmp_pk(tmp); tmp_fr = tmp_fr(tmp); tmp_wd = tmp_wd(tmp);
    if length(tmp_pk) > 3
        tmp_pk = tmp_pk(1:3); tmp_fr = tmp_fr(1:3); tmp_wd = tmp_wd(1:3);
    end
    [~,tmp] = sort(tmp_fr, 'ascend');

    osci.pk.pk{e} = tmp_pk(tmp);
    osci.pk.fr{e} = tmp_fr(tmp);
    osci.pk.wd{e} = tmp_wd(tmp);
    clear tmp
        
    tmp = osci.pk.fr{e} > 4 & osci.pk.fr{e} < 8.5;
    if sum(tmp) == 1
        osci.theta(e) = osci.pk.fr{e}(tmp);
    elseif sum(tmp) > 1
        [~,tmp2] = sort(osci.pk.pk{e}(tmp), 'descend');
        tmp = osci.pk.fr{e}(tmp);
        osci.theta(e) = tmp(tmp2(1));
    end
    clear tmp tmp2
    
    tmp = osci.pk.fr{e} > 8 & osci.pk.fr{e} < 13.5;
    if sum(tmp) == 1
        osci.alpha(e) = osci.pk.fr{e}(tmp);
    elseif sum(tmp) > 1
        [~,tmp2] = sort(osci.pk.pk{e}(tmp), 'descend');
        tmp = osci.pk.fr{e}(tmp);
        osci.alpha(e) = tmp(tmp2(1));
    end
    clear tmp tmp2
    
    tmp = osci.pk.fr{e} > 13;
    if sum(tmp) == 1
        osci.beta(e) = osci.pk.fr{e}(tmp);
    elseif sum(tmp) > 1
        [~,tmp2] = sort(osci.pk.pk{e}(tmp), 'descend');
        tmp = osci.pk.fr{e}(tmp);
        osci.beta(e) = tmp(tmp2(1));
    end
    clear tmp*
end

% save oscillatory peaks, original power spectrum, and fractal component
save(fullfile(savdir, 'osci'), 'osci', 'orig', 'frac');
clear osci orig frac

end
