function [k,velplot] = cellpopsdiff_upwind(j_init,c1_init,c2_init,c,r,t,p1,p2,kappa,mu,alpha,beta,gamma,cmin,rbar)

%%% spatial mesh
R = length(r);
dr = r(2)-r(1);

%%% time mesh
dt = diff(t);
T = length(t);

%%% initialize
k = zeros(T,R);
k(1,:) = c2_init - c1_init;

velplot=zeros(T,R);

for i=2:T
    jj = j_init + (i-2);
    Tp = Tderivative(c(i,:),kappa,cmin,rbar);
    v = zeros(1,R);
    v(2:jj) = 1/mu .* Tp(2:jj) .* (c(i,3:jj+1) - c(i,1:jj-1))/(2*dr);
    v(1) = 0;
    v(jj+1) = 1/mu * Tp(jj+1) * (3/2*c(i,jj+1) - 2*c(i,jj) + 1/2*c(i,jj-1))/dr;
    
    maindiag = 1 + dt(i-1)/dr *v(1:jj+1) - dt(i-1)*( alpha*p1(i,1:jj+1) ...
        - beta*p2(i,1:jj+1) - gamma );
    upperdiag = [2*dt(i-1)/dr * v(2) , zeros(1,jj-1)];
    lowerdiag = -dt(i-1)./dr .* ( r(1:jj)./r(2:jj+1) ) .*v(1:jj);
    
    thetamatrix = diag(maindiag) + diag(upperdiag,1) + diag(lowerdiag,-1);
    
    kextrap = interp1(r(1:jj),k(i-1,1:jj),r(jj+1),'pchip','extrap');
    
    bvector = ( [k(i-1,1:jj) , kextrap] + dt(i-1)*beta.*p2(i,1:jj+1).*c(i,1:jj+1) )';
    
    k_new = ( thetamatrix \ bvector )';

    k(i,:) = [k_new , zeros(1,R-(jj+1))];
    
    velplot(i,:) = v;
%     keyboard
    
%     figure(1)
%     hold on
%     plot(k_new)
%     pause(.5)
    
    
%     r_interp = linspace(0,r(jj),jj+1);
%     c1_interp = [interp1(r(1:jj),c1(i-1,1:jj),r_interp,'pchip') , zeros(1,R-(jj+1))];
    
end