function k_new = cellpops_sum_withgrowthfactors(j,c1_old,c2_old,q1,PO2,...
    dt,r,Pm,kappa,mu,alpha11,alpha12,alpha21,alpha22,gamma1,gamma2,...
    cmin,rbar,ce,cmax,hyaloid)
% k_new = cellpops_sum_withgrowthfactors(j,c1_old,c2_old,q1,PO2,...
%     dt,r,Pm,kappa,mu,alpha11,alpha12,alpha21,alpha22,gamma1,gamma2,...
%     cmin,rbar,ce,cmax,hyaloid)
%
% inputs:
%   j      = node that moving boundary is located at
%   c1_old = APC cell concentration at previous time
%   c2_old = IPA cell concentration at previous time
%   q1     = PDGFA at next time
%   q2     = LIF at next time
%   PO2    = partial pressure of oxygen at next time
%   dt     = time step size
%   r      = spatial mesh
%   {others} = parameters
%
% outputs:
%   k_new = APC+IPA cell concentration at next time
%
% note: cells at time step j are at mesh points 1:j
%       cells at time step j+1 are at mesh points 1:j+1

global whatstep tcurr;

%%% spatial mesh
R = length(r);
dr = r(2)-r(1);
rhalf = (r(1:R-1)+r(2:R))/2;

%%% initialize
k_old = c1_old + c2_old;
khalf = (k_old(1:R-1)+k_old(2:R))/2;
Tp = Tderivative(khalf,kappa,cmin,rbar);
Psi = rhalf.*Tp.*khalf;

%%% extrapolate Psi to get node j+1/2
% Psi(j) = interp1(rhalf(1:j-1),Psi(1:j-1),rhalf(j),'pchip','extrap');
Psi(j) = interp1(rhalf(j-2:j-1),Psi(j-2:j-1),rhalf(j),'pchip','extrap');

% rinterp = linspace(0,r(j),j+1);
% rinterp_half = (rinterp(1:end-1)+rinterp(2:end))/2;
% kinterp = [interp1(r(1:j),k_old(1:j),rinterp_half,'pchip') , zeros(1,R-j-1)];
% Tp = Tderivative(kinterp,kappa,cmin,rbar);
% Psi = rhalf.*Tp.*kinterp;

omega = 1./(mu*r(2:j)) * dt/dr^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% growth function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g = growthterms_sum_withgrowthfactors(c1_old(1:j),c2_old(1:j),q1(1:j),...
    PO2(1:j),Pm,alpha11,alpha12,alpha21,alpha22,gamma1,gamma2,cmax,hyaloid);

%%%%%%%%%%%%%%%%%%%%%%%%% iterate for convergence %%%%%%%%%%%%%%%%%%%%%%%%%
for m=1:5
    %%%%%%%%%%%%%%%%%%%%%%%%%% construct matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%
    upperdiag = [2/mu * dt/dr^2 * khalf(1)*Tp(1) , ...
        omega .* Psi(2:j)];
    maindiag = [1 - 2/mu * dt/dr^2 * khalf(1)*Tp(1) , ...
        1 - omega .* ( Psi(2:j) + Psi(1:j-1) ) , ...
        1];
    % upperdiag = [2/mu * dt/dr^2 * kinterp(1)*Tp(1) , ...
    %     omega .* Psi(2:j)];
    % maindiag = [1 - 2/mu * dt/dr^2 * kinterp(1)*Tp(1) , ...
    %     1 - omega .* ( Psi(2:j) + Psi(1:j-1) ) , ...
    %     1];
    lowerdiag = [omega .* Psi(1:j-1) , ...
        0];

    thetamatrix = diag(maindiag) + diag(upperdiag,1) + diag(lowerdiag,-1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%% right hand side %%%%%%%%%%%%%%%%%%%%%%%%%%%
    bvector = [k_old(1:j)' + dt*g' ; ce];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%% solve system %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    k_new = ( thetamatrix \ bvector )';

    k_new = [k_new(1:j+1) , zeros(1,R-(j+1))];

    %%%%%%%%%%%%%%%%%%%%%% update for next iteration %%%%%%%%%%%%%%%%%%%%%%
    khalf = (k_new(1:R-1)+k_new(2:R))/2;
    Tp = Tderivative(khalf,kappa,cmin,rbar);
    Psi = rhalf.*Tp.*khalf;

%     if tcurr>3.22*24
%         hold on
%         if dr>0.01
%             plot(r,k_new,'-o')
%         else
%             plot(r,k_new)
%         end
%         hold off
%     end

end

% if tcurr>3.22*24
%     disp(whatstep)
%     keyboard
% end

% -1/(mu*r(j)*dr^2) * ...
% ( rhalf(j)*(k_new(j)+k_new(j+1))/2*Tderivative((k_new(j)+k_new(j+1))/2,kappa,cmin,rbar)*(k_new(j+1)-k_new(j)) ...
% - rhalf(j-1)*(k_new(j-1)+k_new(j))/2*Tderivative((k_new(j-1)+k_new(j))/2,kappa,cmin,rbar)*(k_new(j)-k_new(j-1)) ) ...
% +g(j)
% 
% (k_new(j)-k_old(j))/dt
% 
% keyboard