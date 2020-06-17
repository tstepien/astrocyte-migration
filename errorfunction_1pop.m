function [err_tot,err_time,err_rad] = errorfunction_1pop(t,mvgbdy)
% [err_tot,err_time,err_rad] = errorfunction_1pop(t,mvgbdy)
%
% This is the error function for comparing the experimental data with
% simulations
%
% inputs:
%   t      = time vector
%   mvgbdy = vector with location of moving boundary over time
%
% outputs:
%   err_tot  = total error
%   err_time = error from time end point
%   err_rad  = error from astrocyte radius

%%% penalties
if t(end)/24>8  || ~isreal(t(end)) %|| t(end)/24 <6
    err_tot = 10^4; %NaN;
    err_time = 10^4; %NaN;
    err_rad = 10^4; %NaN;
    return
end

%%% fixed parameters
parameters_fixed

%%% APC radius (mm) for E15-E16, E18-E22/P0
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
    ind(i) = find(abs((t/24-(jj-1)))==min(abs(t/24-(jj-1))),1,'first');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% errors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% time error
err_time = abs(7 - t(end)/24)/7;

%%% radius error
err_rad = sum( abs(rad_APC - mvgbdy(ind)) ./ rad_APC );

%%% total error
err_tot = err_time + err_rad;