function err_rad = errorfunction_1pop(t,mvgbdy)
% err_rad = errorfunction_1pop(t,mvgbdy)
%
% This is the error function for comparing the experimental data with
% simulations
%
% inputs:
%   t      = time vector
%   mvgbdy = vector with location of moving boundary over time
%
% outputs:
%   err_rad  = error from astrocyte radius

if t(end)/24 <6 || t(end)/24>8
    err_rad = NaN;
    return
end

%%% threshold for small density
parameters_fixed

%%% APC, IPA, and retina radius (mm) for E15-E16, E18-E22/P0
rad_APC = [0.17; 0.33; 0.5; 0.67; 1.67; 2.17; 2.67];

%%% times: day 0, 1, 2, 3, 4, 5, 6, 7
%%% but note that we don't have data for day 2
dayswithdata = [1:2 4:8];
numdays = length(dayswithdata);

%%% pre-allocate arrays
ind = zeros(numdays,1);

%%% calculate values on each day of data
for i=1:numdays
    jj = dayswithdata(i);
    ind(i) = find(abs((t/24-(jj-1)))==min(abs(t/24-(jj-1))));
end

%%% total error from astrocyte radius
err_rad = sum( abs(rad_APC - mvgbdy(ind)) ./ rad_APC );