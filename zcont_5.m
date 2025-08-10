function z_mtx = zcont_5(data_mtx)
% ZBASELINE_200 - subfunction to z-score continous data on itself, randomly
% drawing 5 data points per permutation.
%
% Inputs:
% data = data matrix with 2D (i.e., channels x frequencies), to be z-scored
% 
% Outputs:
% z_mtx = z-scored data matrix that is otherwise identical to the input 
%   data matrix
%
% Copyright (c) 2023
% EL Johnson, PhD

disp(' '); disp('Performing statistical bootstrapping...'); disp(' ');

npermutes = 1000;
nrand = 5;

z_mtx = nan(size(data_mtx));

ntimes = size(z_mtx,2);

for e = 1:size(z_mtx,1)
    tmp = data_mtx(e,:);

    tmp_dist = zeros(1,npermutes);
    for z = 1:npermutes
        tmp_rand = randsample(tmp, nrand);
        tmp_dist(z) = mean(tmp_rand); 
    end
    
    tmp_mean = repmat(mean(tmp_dist), [1 ntimes]);
    tmp_sd = repmat(std(tmp_dist), [1 ntimes]);

    z_mtx(e,:) = (data_mtx(e,:)-tmp_mean)./tmp_sd;
    clear tmp*
end

clearvars -except z_mtx

end