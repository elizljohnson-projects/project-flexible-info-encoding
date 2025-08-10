function tfrz = zbaseline_200(data, baseline)
% ZBASELINE_200 - subfunction to z-score time-series data on the baseline 
% using statistical bootstrapping, randomly drawing 200 data points per
% permulation.
%
% Inputs:
% data = data structure in standard FieldTrip format containing either 3D 
%   data matrix in 'trial' field (i.e., trials x channels x time) or 4D 
%   data matrix in 'powspctrm' field (i.e., trials x channels x frequencies
%   x time) (e.g., data_stim), to be baseline-normalized
% baseline = baseline data matrix with the same dimensions as the data
%   matrix (e.g., base_mtx)
% 
% Outputs:
% tfrx = baseline-normalized data structure that is otherwise identical to 
%   the input data structure
%
% Copyright (c) 2023
% EL Johnson, PhD

disp(' '); disp('Performing statistical bootstrapping...'); disp(' ');

npermutes = 1000;
nrand = 200;

tfrz = data;

if length(size(baseline)) == 3
    tfrz.trial = zeros(size(tfrz.trial));

    ntrials = size(tfrz.trial, 1);
    ntimes = size(tfrz.trial, 3);

    for e = 1:size(baseline,2)
        base = baseline(:,e,:);
        base = base(:);

        basedist = zeros(1,npermutes);
        for z = 1:npermutes
            brand = randsample(base, nrand);
            basedist(z) = mean(brand); 
        end

        bmean = repmat(mean(basedist), [ntrials 1 ntimes]);
        bsd = repmat(std(basedist), [ntrials 1 ntimes]);

        tfrz.trial(:,e,:) = (data.trial(:,e,:)-bmean)./bsd;

        clear base basedist brand bmean bsd
    end
    
elseif length(size(baseline)) == 4
    tfrz.powspctrm = zeros(size(tfrz.powspctrm));

    ntrials = size(tfrz.powspctrm, 1);
    ntimes = size(tfrz.powspctrm, 4);

    for e = 1:size(baseline,2)
        for f = 1:size(baseline,3)
            base = baseline(:,e,f,:);
            base = base(:);

            basedist = zeros(1,npermutes);
            for z = 1:npermutes
                brand = randsample(base, nrand);
                basedist(z) = mean(brand); 
            end

            bmean = repmat(mean(basedist), [ntrials 1 1 ntimes]);
            bsd = repmat(std(basedist), [ntrials 1 1 ntimes]);

            tfrz.powspctrm(:,e,f,:) = (data.powspctrm(:,e,f,:)-bmean)./bsd;

            clear base basedist brand bmean bsd
        end
    end

clearvars -except tfrz

end