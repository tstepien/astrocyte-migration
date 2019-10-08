function [err_rad,err_dens,err_time,err_tot] = errorfunction(t,r,mvgbdy,c1,c2)
% [err_rad,err_dens,err_time,err_tot] = errorfunction(t,mvgbdy,c1,c2)
%
% This is the error function for comparing the experimental data with
% simulations
%
% inputs:
%   t      = time vector
%   r      = spatial grid vector
%   mvgbdy = vector with location of moving boundary over time
%   c1     = density of APCs
%   c2     = density of IPAs
%
% outputs:
%   err_rad  = error from astrocyte radius
%   err_dens = error from astrocyte density
%   err_time = error from ending time compared to 7 days
%   err_tot  = total error = (err_rad + err_dens + err_time)

%%% APC and IPA radius (mm) for E15-E16, E18-E22/P0
rad_APC = [0.33; 0.33; 0.5; 0.67; 1.67; 2.17; 2.67];
rad_IPA = [0; 0.04; 0.08; 0.33; 1; 1.5; 2];

%%% times: day 0, 1, 2, 3, 4, 5, 6, 7
%%% but note that we don't have data for day 2
dayswithdata = [1:2 4:8];
numdays = length(dayswithdata);

%%% pre-allocate arrays
ind = zeros(numdays,1);
est_APC = zeros(numdays,length(r));
est_IPA = zeros(numdays,length(r));
act_APC = zeros(numdays,length(r));
act_IPA = zeros(numdays,length(r));
numnodes = zeros(numdays,1);

%%% calculate values on each day of data
for i=1:numdays
    jj = dayswithdata(i);
    ind(i) = find(abs((t/24-(jj-1)))==min(abs(t/24-(jj-1))));
    
    est_APC(i,:) = (c1(ind(i),:)>0);
    est_IPA(i,:) = (c2(ind(i),:)>0);
    
    act_APC(i,:) = (r<=rad_APC(i));
    act_IPA(i,:) = (r<=rad_IPA(i));
    
    numnodes(i) = sum(r<=mvgbdy(ind(i)));
end

%%% total error from astrocyte radius
err_rad = sum( abs(rad_APC - mvgbdy(ind)) ./ rad_APC );

%%% total error from density
diff_APC = sum( abs(act_APC - est_APC) , 2) ./ numnodes;
diff_IPA = sum( abs(act_IPA - est_IPA) , 2) ./ numnodes;

err_dens = sum( diff_APC + diff_IPA );

%%% total error from time
err_time = abs(7 - t(end)/24)/7;

%%% total error
err_tot = err_rad + err_dens + err_time;