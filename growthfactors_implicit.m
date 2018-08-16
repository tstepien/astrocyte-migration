function [p1_new,p2_new] = growthfactors_implicit(p1_old,p2_old,dt,r,...
    D1,D2,eta1,eta2)
% [p1_new,p2_new] = growthfactors_implicit(p1_old,p2_old,dt,r,D1,D2,eta1,eta2)
%
% inputs:
%   p1_old = PDGFA growth factor concentration at previous time
%   p2_old = LIF growth factor concentration at previous time
%   dt     = time step size
%   r      = spatial mesh
%   {others} = parameters
%
% outputs:
%   p1_new = PDGFA growth factor concentration at next time
%   p2_new = LIF growth factor concentration at next time

%%% spatial mesh
R = length(r);
dr = r(2)-r(1);

%%% Neumann boundary conditions at end of domain (r=rmax)
%%% (partial p/partial t) = constant
p1BC = 0;
p2BC = 0;

%%% iterate based off of CFL condition
if dt >= dr^2/(2*D1)
% %     newdt = (dr^2/(2*D1)) *(7/8);
% %     t = 0:newdt:dt;
% %     if t(end)<dt
% %         t = [t, dt];
% %     end
% %     num_iter = length(t)-1;
% %     tnew = linspace(0,dt,num_iter+1);
% %     dt = tnew(2)-tnew(1);

    %%% /\ above iteration will result in waaaay too many iterations
    %%% \/ below iteration controls that we only add on 10 iterations to
    %%%    attempt to control blow up (which shouldn't be happening...)
    num_iter = 5;
    dt = dt/num_iter;
else
    num_iter = 1;
end

for i=1:num_iter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% p1 - PDGFA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta1_1 = -D1 * dt/dr^2 * (1 + dr./(2*r(2:R-1)));
theta2_1 = 1 + 2*D1*dt/dr^2*ones(1,R-1) - dt*eta1;
theta3_1 = -D1 * dt/dr^2 * (1 - dr./(2*r(2:R-1)));

theta4_1 = -4*D1 * dt/dr^2;
theta5_1 = 1 + 4*D1 * dt/dr^2 - dt*eta1;
theta6_1 = -2*D1 * dt/dr^2;
theta7_1 = 2*D1 * dt/dr^2 * p1BC*dr * ( 1 + dr/(2*r(R)) );

maindiag1 = [theta5_1 , theta2_1];
upperdiag1 = [theta4_1 , theta1_1];
lowerdiag1 = [theta3_1 , theta6_1];

thetamatrix1 = diag(maindiag1) + diag(upperdiag1,1) + diag(lowerdiag1,-1);

bvector1 = p1_old' + [zeros(R-1,1) ; theta7_1];

p1_new = ( thetamatrix1 \ bvector1 )';
p1_old = p1_new;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% p2 - LIF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta1_2 = -D2 * dt/dr^2 * (1 + dr./(2*r(2:R-1)));
theta2_2 = 1 + 2*D2*dt/dr^2*ones(1,R-1) - dt*eta2;
theta3_2 = -D2 * dt/dr^2 * (1 - dr./(2*r(2:R-1)));
theta4_2 = -4*D2 * dt/dr^2;
theta5_2 = 1 + 4*D2 * dt/dr^2 - dt*eta2;
theta6_2 = 0;%-2*D2 * dt/dr^2;
theta7_2 = 0;%2*D2 * dt/dr^2 * p2BC*dr * ( 1 + dr/(2*r(R)) );

maindiag2 = [theta5_2 , theta2_2];
upperdiag2 = [theta4_2 , theta1_2];
lowerdiag2 = [theta3_2 , theta6_2];

thetamatrix2 = diag(maindiag2) + diag(upperdiag2,1) + diag(lowerdiag2,-1);
thetamatrix2(end,end) = 1;

bvector2 = p2_old' + [zeros(R-1,1) ; theta7_2];

p2_new = ( thetamatrix2 \ bvector2 )';
p2_old = p2_new;

end