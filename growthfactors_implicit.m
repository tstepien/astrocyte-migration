function [q1_new,q2_new] = growthfactors_implicit(q1_old,q2_old,dt,tcurr,...
    r,dr,R,thickness_RGC,radius_endo,maxRGCthick,thickness_ret,D1,D2,xi1,...
    xi2,gamma3,gamma4)
% [q1_new,q2_new] = growthfactors_implicit(q1_old,q2_old,dt,tcurr,...
%     r,dr,R,thickness_RGC,radius_endo,maxRGCthick,thickness_ret,D1,D2,xi1,...
%     xi2,gamma3,gamma4)
%
% Uses Dirichlet boundary condition q1=q2=0 at Rmax
%
% inputs:
%   q1_old = PDGFA growth factor concentration at previous time
%   q2_old = LIF growth factor concentration at previous time
%   dt     = time step size
%   tcurr  = current time
%   r      = spatial mesh
%   dr     = spatial mesh size
%   R      = length of spatial mesh
%   {others} = parameters
%
% outputs:
%   q1_new = PDGFA growth factor concentration at next time
%   q2_new = LIF growth factor concentration at next time


tday = (tcurr+dt)/24;
if tday<3
    timeramp = 0;
else
    timeramp = 1/4*(tday-3);
end

if radius_endo==0
    spaceind = (r<radius_endo)';
else
    spaceind = (r<=radius_endo)';
end

%%% spatial mesh
% nodesretina = sum(thickness_ret>0); %%% PDGFA and LIF can spread through the
%                                 %%% current extent of the retina
% Rorig = length(r);
% r = r(1:nodesretina); %%% restrict domain
% q1_old = q1_old(1:nodesretina); %%% restrict domain
% q2_old = q2_old(1:nodesretina);

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
theta2_1 = 1 + 2*D1*dt/dr^2*ones(1,R-1) + dt*gamma3;%*thickness_RGC(2:R)/maxRGCthick;
theta1_1 = -D1 * dt/dr^2 * (1 - dr./(2*r(2:R-1)));

% origin
theta5_1 = -4*D1 * dt/dr^2;
theta4_1 = 1 + 4*D1 * dt/dr^2 + dt*gamma3;%*thickness_RGC(1)/maxRGCthick;

% rmax
theta6_1 = 0;%-2*D1 * dt/dr^2;
theta7_1 = 0;%2*D1 * dt/dr^2 * p1BC*dr * ( 1 + dr/(2*r(R)) );

maindiag1 = [theta4_1 , theta2_1];
upperdiag1 = [theta5_1 , theta3_1];
lowerdiag1 = [theta1_1 , theta6_1];

thetamatrix1 = diag(maindiag1) + diag(upperdiag1,1) + diag(lowerdiag1,-1);
thetamatrix1(end,end) = 1;

bvector1 = q1_old' + [zeros(R-1,1) ; theta7_1] ...
    + dt*xi1.*thickness_RGC'/maxRGCthick * timeramp.*(thickness_ret>0)';

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

bvector2 = q2_old' + [zeros(R-1,1) ; theta7_2] + dt*xi2.*spaceind;

q2_new = ( thetamatrix2 \ bvector2 )';
% % % q2_old = q2_new;

% % % end

%%% resize
% q1_new = [q1_new , zeros(1,Rorig-nodesretina)];
% q2_new = [q2_new , zeros(1,Rorig-nodesretina)];