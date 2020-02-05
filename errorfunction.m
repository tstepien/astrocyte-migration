function [err_rad,err_dens,err_tot] = errorfunction(t,r,mvgbdy,c1,c2)
% [err_rad,err_dens,err_tot] = errorfunction(t,mvgbdy,c1,c2)
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
%   err_tot  = total error = (err_rad + err_dens + err_time)

if t(end)/24 <6 || t(end)/24>8
    err_rad = NaN;
    err_dens = NaN;
    err_tot = NaN;
    return
end

%%% threshold for small density
parameters_fixed
thd_low = 50;
thd_high = cmin;

%%% APC, IPA, and retina radius (mm) for E15-E16, E18-E22/P0
rad_APC = [0.33; 0.33; 0.5; 0.67; 1.67; 2.17; 2.67];
rad_IPA = [0; 0.04; 0.08; 0.33; 1; 1.5; 2];
rad_ret = [1.17; 1.5; 2.17; 2.5; 2.83; 3.83; 4];

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
dens_annulus = zeros(numdays,length(r));
dens_disc = zeros(numdays,length(r));
% numnodes = zeros(numdays,1);
nodes_APC = zeros(numdays,length(r));
nodes_IPA = zeros(numdays,length(r));
numnodes_APC = zeros(numdays,1);
numnodes_IPA = zeros(numdays,1);

%%% calculate values on each day of data
for i=1:numdays
    jj = dayswithdata(i);
    ind(i) = find(abs((t/24-(jj-1)))==min(abs(t/24-(jj-1))));
    
    est_APC(i,:) = (c1(ind(i),:)>0);
    est_IPA(i,:) = (c2(ind(i),:)>0);
    
%     numnodes(i) = sum(r<=mvgbdy(ind(i)));
    nodes_APC(i,:) = (r<=rad_APC(i) & r>=rad_IPA(i)); %% annulus
    nodes_IPA(i,:) = (r<rad_IPA(i)); %% disc
    
    numnodes_APC(i) = sum(nodes_APC(i,:));
    numnodes_IPA(i) = sum(nodes_IPA(i,:));
    
    %%% incorrect density relationship for APCs and IPAs on (APC annulus)
    %%% and (IPA disc)
    dens_annulus(i,:) = (1 - ( c2(ind(i),:) < cmin & ...
        c1(ind(i),:)>cmin ) ) .* nodes_APC(i,:);
    dens_disc(i,:) = (1 - ( c1(ind(i),:) < cmin & ...
        c2(ind(i),:)>cmin ) ) .* nodes_IPA(i,:);
    
%     act_APC(i,:) = (r<=rad_APC(i) & r>rad_IPA(i)); %% annulus
%     act_IPA(i,:) = (r<=rad_IPA(i)); %% disc
end

%%% total error from astrocyte radius
err_rad = sum( abs(rad_APC - mvgbdy(ind)) ./ rad_APC );

%%% total error from density
err_APC = sum( dens_annulus , 2) ./ numnodes_APC;
err_IPA = sum( dens_disc , 2) ./ numnodes_IPA;
err_IPA(1) = 0; %%% initial time point has no IPAs by initial condition, so
                %%% numnodes_IPA=0 and dividing by zero results in NaN
% diff_APC = sum( abs(act_APC - est_APC) , 2) ./ numnodes_APC;
% diff_IPA = sum( abs(act_IPA - est_IPA) , 2) ./ numnodes_IPA;

err_dens = sum( err_APC + err_IPA );

%%% total error
err_tot = err_rad + err_dens;

if sum(sum(c1<0))>0 || sum(sum(c2<0))>0
    err_tot = NaN;
end

% figure
% plot(r,dens_disc)
% title('(incorrect disc) where c2<c1')
% figure
% plot(r,dens_annulus)
% title('(incorrect annulus) where c1<c2')
% figure
% plot(r,nodes_IPA)
% title('nodes where IPAs are (disc)')
% figure
% plot(r,nodes_APC)
% title('nodes where APCs are (annulus)')
% 
% keyboard