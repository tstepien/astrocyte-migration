function [q1_new,q2_new] = growthfactors_implicit(q1_old,q2_old,dt,tcurr,...
    r,D1,D2,xi1,xi2,gamma3,gamma4,thickness,width_retina)
% [q1_new,q2_new] = growthfactors_implicit(q1_old,q2_old,dt,tcurr,...
%     r,D1,D2,xi1,xi2,gamma3,gamma4,thickness,width_retina)
%
% inputs:
%   q1_old = PDGFA growth factor concentration at previous time
%   q2_old = LIF growth factor concentration at previous time
%   dt     = time step size
%   tcurr  = current time
%   r      = spatial mesh
%   {others} = parameters
%
% outputs:
%   q1_new = PDGFA growth factor concentration at next time
%   q2_new = LIF growth factor concentration at next time

%%% CURRENTLY USES DIRICHLET=0 BC AT RMAX

%%% convert current time from hours to days
tday = (tcurr+dt)/24;
timeind = (tday>=3); %day 3 = E18

%%% spatial mesh
% nodesretina = sum(thickness>0); %%% PDGFA and LIF can spread through the
%                                 %%% current extent of the retina
% Rorig = length(r);
% r = r(1:nodesretina); %%% restrict domain
R = length(r);
dr = r(2)-r(1);

% q1_old = q1_old(1:nodesretina); %%% restrict domain
% q2_old = q2_old(1:nodesretina);

%%% thickness of retinal ganglion cell layer (from Braekevelt and Hollenburg)
%%% in microns, converted to mm
thickness_posterior = max(-1.95*tday^2 + 14.84*tday + 9.01 , 0) * 0.001;
thickness_peripheral = max(-3.66*tday^2 + 26.83*tday - 14.7 , 0) * 0.001;

thickness_RGC = max( (thickness_peripheral-thickness_posterior)/width_retina^2 ...
    * r.^2 + thickness_posterior, 0) .* (r<=width_retina);
maxthick = 46 * 0.001; % maximum thickness of RGC layer: 46 micon converted to mm

[sum(thickness>0) sum(thickness_RGC>0)];

%%% radius of endothelial cells (in microns, converted to mm)
radius_endo = max(425*tday - 1675 , 0) * 0.001;

%%% Neumann boundary conditions at end of domain (r=rmax)
%%% (partial p/partial t) = constant
% p1BC = 0;

% % % %%% iterate based off of CFL condition
% % % if dt >= dr^2/(2*D1)
% % % % %     newdt = (dr^2/(2*D1)) *(7/8);
% % % % %     t = 0:newdt:dt;
% % % % %     if t(end)<dt
% % % % %         t = [t, dt];
% % % % %     end
% % % % %     num_iter = length(t)-1;
% % % % %     tnew = linspace(0,dt,num_iter+1);
% % % % %     dt = tnew(2)-tnew(1);
% % % 
% % %     %%% /\ above iteration will result in waaaay too many iterations
% % %     %%% \/ below iteration controls that we only add on 10 iterations to
% % %     %%%    attempt to control blow up (which shouldn't be happening...)
% % %     num_iter = 5;
% % %     dt = dt/num_iter;
% % % else
% % %     num_iter = 1;
% % % end

% % % for i=1:num_iter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% q1 - PDGFA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta3_1 = -D1 * dt/dr^2 * (1 + dr./(2*r(2:R-1)));
theta2_1 = 1 + 2*D1*dt/dr^2*ones(1,R-1) + dt*gamma3;%*thickness_RGC(2:R)/maxthick;
theta1_1 = -D1 * dt/dr^2 * (1 - dr./(2*r(2:R-1)));

% origin
theta5_1 = -4*D1 * dt/dr^2;
theta4_1 = 1 + 4*D1 * dt/dr^2 + dt*gamma3;%*thickness_RGC(1)/maxthick;

% rmax
theta6_1 = 0;%-2*D1 * dt/dr^2;
theta7_1 = 0;%2*D1 * dt/dr^2 * p1BC*dr * ( 1 + dr/(2*r(R)) );

maindiag1 = [theta4_1 , theta2_1];
upperdiag1 = [theta5_1 , theta3_1];
lowerdiag1 = [theta1_1 , theta6_1];

thetamatrix1 = diag(maindiag1) + diag(upperdiag1,1) + diag(lowerdiag1,-1);
thetamatrix1(end,end) = 1;

bvector1 = q1_old' + [zeros(R-1,1) ; theta7_1] ...
    + dt*xi1.*thickness_RGC'/maxthick * timeind.*(thickness>0)';

q1_new = ( thetamatrix1 \ bvector1 )';
% % % q1_old = q1_new;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% q2 - LIF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta3_2 = -D2 * dt/dr^2 * (1 + dr./(2*r(2:R-1)));
theta2_2 = 1 + 2*D2*dt/dr^2*ones(1,R-1) + dt*gamma4;
theta1_2 = -D2 * dt/dr^2 * (1 - dr./(2*r(2:R-1)));

% origin
theta5_2 = -4*D2 * dt/dr^2;
theta4_2 = 1 + 4*D2 * dt/dr^2 + dt*gamma4;

% rmax
theta6_2 = 0;
theta7_2 = 0;

maindiag2 = [theta4_2 , theta2_2];
upperdiag2 = [theta5_2 , theta3_2];
lowerdiag2 = [theta1_2 , theta6_2];

thetamatrix2 = diag(maindiag2) + diag(upperdiag2,1) + diag(lowerdiag2,-1);
thetamatrix2(end,end) = 1;

bvector2 = q2_old' + [zeros(R-1,1) ; theta7_2] + dt*xi2.*(r<=radius_endo)';

q2_new = ( thetamatrix2 \ bvector2 )';
% % % q2_old = q2_new;

% % % end

%%% resize
% q1_new = [q1_new , zeros(1,Rorig-nodesretina)];
% q2_new = [q2_new , zeros(1,Rorig-nodesretina)];

% keyboard