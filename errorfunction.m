function [err_tot,err_time,err_rad,err_dens] = errorfunction(t,r,mvgbdy,c1,c2)
% [err_tot,err_time,err_rad,err_dens] = errorfunction(t,mvgbdy,c1,c2)
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

%%% penalties
if t(end)/24>8  || ~isreal(t(end)) || sum(c1(:)<0)>0 ...
        || sum(c2(:)<0)>0 %|| t(end)/24 <6
    err_tot = 10^4; %NaN;
    err_time = 10^4; %NaN;
    err_rad = 10^4; %NaN;
    err_dens = 10^4; %NaN;
    return
end

%%% fixed parameters
parameters_fixed

%%% APC and IPA radius (mm) for E15-E16, E18-E22/P0
rad_APC = [0.17; 0.33; 0.5; 0.67; 1.67; 2.17; 2.67];
rad_IPA = [0; 0.04; 0.08; 0.33; 1; 1.5; 2];

%%% times: day 0, 1, 2, 3, 4, 5, 6, 7
%%% but note that we don't have data for day 2
dayswithdata = [1:2 4:8];
numdays = length(dayswithdata);

%%% pre-allocate arrays
ind = zeros(numdays,1);
comp_APC = zeros(numdays,length(r));
comp_IPA = zeros(numdays,length(r));
nodes_APC = zeros(numdays,length(r));
nodes_IPA = zeros(numdays,length(r));
numnodes_APC = zeros(numdays,1);
numnodes_IPA = zeros(numdays,1);
dens_annulus = zeros(numdays,length(r));
dens_disc = zeros(numdays,length(r));

%%% calculate values on each day of data
for i=1:numdays
    %%% radius
    jj = dayswithdata(i);
    ind(i) = find(abs((t/24-(jj-1)))==min(abs(t/24-(jj-1))),1,'first');
    
    %%% density
    comp_APC(i,:) = (c1(ind(i),:)>0);
    comp_IPA(i,:) = (c2(ind(i),:)>0);
    
    nodes_APC(i,:) = (r<=rad_APC(i) & r>rad_IPA(i)); %% annulus
    nodes_IPA(i,:) = (r<=rad_IPA(i)); %% disc
    
    numnodes_APC(i) = sum(nodes_APC(i,:));
    numnodes_IPA(i) = sum(nodes_IPA(i,:));
    
    %%% incorrect density relationship for APCs and IPAs on
    %%% (APC annulus) and (IPA disc)
    
    %%% where c1<cmin or c2>cmin on the APC annulus (incorrect relationship)
    dens_annulus(i,:) = ( c2(ind(i),:)>cmin | c1(ind(i),:)<cmin ) .* nodes_APC(i,:);
    %%% where c1>cmin or c2<cmin on the IPA disc (incorrect relationship)
    dens_disc(i,:) = ( c1(ind(i),:)>cmin | c2(ind(i),:)<cmin ) .* nodes_IPA(i,:);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% errors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% time error
err_time = abs(7 - t(end)/24)/7;

%%% radius error
err_rad = sum( abs(rad_APC - mvgbdy(ind)) ./ rad_APC );

%%% density error
err_APC = sum( dens_annulus , 2) ./ numnodes_APC;
err_IPA = sum( dens_disc , 2) ./ numnodes_IPA;
err_IPA(1) = 0; %%% initial time point has no IPAs by initial condition, so
                %%% numnodes_IPA=0 and dividing by zero results in NaN
err_dens = sum( err_APC + err_IPA );

%%% total error
err_tot = err_time + err_rad + err_dens;

%%% figure
% figure
% tiledlayout(2,2)
% 
% nexttile
% plot(r,dens_annulus)
% title('where c1<cmin or c2>cmin on the APC annulus (incorrect relationship)')
% 
% nexttile
% plot(r,dens_disc)
% title('where c1>cmin or c2<cmin on the IPA disc (incorrect relationship)')
% 
% nexttile
% plot(r,nodes_APC)
% title('nodes where APCs are (annulus)')
% 
% nexttile
% plot(r,nodes_IPA)
% title('nodes where IPAs are (disc)')