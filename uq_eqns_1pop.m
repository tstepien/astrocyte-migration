function Y = uq_eqns_1pop(X)
% Y = uq_eqns_1pop(X)
%
% This function returns the value of the error of astrocyte migration for a
% single population of cells, described by 7 variables
%
% inputs:
%   X = [mu, alpha11, alpha12, gamma1, Te, P_hy, r_hy]
%
% output:
%   Y = total error

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
p.alpha21 = 0;
p.alpha22 = 0;
p.beta = 0;
p.beta_hat = 0;
p.gamma2 = 0;

%%% mesh parameters
m.dr = 0.01;
m.rmax = 5;
m.tmax = 7*24;

%%% plots on or off
plotsonoff = 'off';

%%%%%%%%%%%%%%%%%%% solve equation and calculate error %%%%%%%%%%%%%%%%%%%%
%%% initialize
N = size(X,1);
Y = zeros(N,1);

for i=1:N
    p.mu = mu(i);
    p.alpha11 = alpha11(i);
    p.alpha12 = alpha12(i);
    p.gamma1 = gamma1(i);
    p.Te = Te(i);
    p.P_hy = P_hy(i);
    p.r_hy = r_hy(i);
    
    [t,~,~,~,~,~,mvgbdy,~,~] = eqnsolver(p,m,plotsonoff);
    
    [Y(i),~,~] = errorfunction_1pop(t,mvgbdy);
end