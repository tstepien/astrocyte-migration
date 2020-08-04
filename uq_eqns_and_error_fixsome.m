function Y = uq_eqns_and_error_fixsome(X)
% Y = uq_eqns_and_error_fixsome(X)
%
% This function returns the computed values of the moving boundary location
% astrocyte migration for a single population of cells and the final
% simulated time (days)
%
% inputs:
%   X = [mu, alpha11, alpha12, alpha22, beta1, beta3, gamma2, P_hy, r_hy]
%
% output:
%   Y = total error

%%%%%%%%%%%%%%%%%%%%%%%%%% parameters to examine %%%%%%%%%%%%%%%%%%%%%%%%%%
%%% APCs and IPAs
mu = X(:,1);
alpha11 = X(:,2);
alpha12 = X(:,3);
alpha22 = X(:,4);
beta1 = X(:,5);
beta3 = X(:,6);
gamma2 = X(:,7);

%%% hyaloid artery
P_hy = X(:,8);
r_hy = X(:,9);

%%% fixed parameters
alpha21 = 0;
beta2 = 0;
gamma1 = 0;
Te = 0.0035;

%%%%%%%%%%%%%%%%%%%%%%%%%%%% fixed parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% mesh parameters
m.dr = 0.01;
m.rmax = 5;
m.tmax = 7*24;

%%%%%%%%%%%%%%%%%%% solve equation and calculate error %%%%%%%%%%%%%%%%%%%%
%%% initialize
N = size(X,1);
Y = zeros(N,1);

for i=1:N
    [t,r,c1,c2,~,~,mvgbdy,~,~] = eqnsolver(mu(i),alpha11(i),alpha12(i),...
        alpha21,alpha22(i),beta1(i),beta2,beta3(i),gamma1,...
        gamma2(i),Te,P_hy(i),r_hy(i),m);
    
    [Y(i),~,~,~] = errorfunction(t,r,mvgbdy,c1,c2);
end