function Y = uq_eqns_1pop(X)
% Y = uq_eqns_1pop(X)
%
% This function returns the computed values of the moving boundary location
% astrocyte migration for a single population of cells and the final
% simulated time (days)
%
% inputs:
%   X = [mu, alpha11, alpha12, gamma1, Te, P_hy, r_hy]
%
% output:
%   Y = moving boundary location at data point times + ending time (days)

%%% times: day 0, 1, 2, 3, 4, 5, 6, 7
%%% but note that we don't have data for day 2
dayswithdata = [1:2 4:8];
numdays = length(dayswithdata);

%%%%%%%%%%%%%%%%%%%%%%%%%% parameters to examine %%%%%%%%%%%%%%%%%%%%%%%%%%
%%% only APCs
mu = X(:,1);
alpha11 = X(:,2);
alpha12 = X(:,3);
gamma1 = X(:,4);

%%% moving boundary tension
Te = X(:,5);

%%% hyaloid artery
P_hy = X(:,6);
r_hy = X(:,7);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% fixed parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% set all IPA stuff = 0
alpha21 = 0;
alpha22 = 0;
beta = 0;
beta_hat = 0;
gamma2 = 0;

%%% mesh parameters
m.dr = 0.01;
m.rmax = 5;
m.tmax = 7*24;

%%%%%%%%%%%%%%%%%%% solve equation and calculate error %%%%%%%%%%%%%%%%%%%%
%%% initialize
N = size(X,1);
Y = zeros(N,numdays+1);

parfor i=1:N
    [t,~,~,~,~,~,mvgbdy,~,~] = eqnsolver(mu(i),alpha11(i),alpha12(i),...
        alpha21,alpha22,beta,beta_hat,gamma1(i),gamma2,Te(i),P_hy(i),r_hy(i),m);
    
    ind = zeros(numdays,1);
    for k=1:numdays
        jj = dayswithdata(k);
        ind(k) = find(abs((t/24-(jj-1)))==min(abs(t/24-(jj-1))),1,'first');
    end
    
    Y(i,:) = [mvgbdy(ind)' t(end)/24];
end