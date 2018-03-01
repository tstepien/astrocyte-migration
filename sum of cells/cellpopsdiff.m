function k = cellpopsdiff(j_init,c1_init,c2_init,c,r,t,p1,p2,kappa,mu,alpha,beta,gamma,cmin,rbar)

%%% spatial mesh
R = length(r);
dr = r(2)-r(1);

%%% time mesh
dt = diff(t);
T = length(t);

%%% initialize
k = zeros(T,R);
k(1,:) = c2_init - c1_init;

for i=2:T
    jj = j_init + (i-2);
    chalf = (c(i,1:R-1)+c(i,2:R))/2;
    Tp = Tderivative(c(i,:),kappa,cmin,rbar);
    v = zeros(1,R);
    v(2:R-1) = 1/mu .* Tp(2:end-1) .* diff(chalf)/dr;
    v(1) = 0;
    v(jj+1) = interp1(r(1:jj),v(1:jj),r(jj+1),'linear','extrap');
    
    theta1 = dt(i-1)./(2*dr*r(2:jj)) .* r(3:jj+1) .* v(3:jj+1);
    theta2 = 1 - dt(i-1)*( alpha*p1(i,1:jj) - beta*p2(i,1:jj) - gamma );
    theta3 = -dt(i-1)./(2*dr*r(2:jj)) .* r(1:jj-1) .* v(1:jj-1);
    
    theta4 = dt(i-1)/dr * v(2);
    
    theta5 = 1 + 3*dt(i-1)/(2*dr) * v(jj+1) - dt(i-1)*( alpha*p1(i,jj+1) ...
        - beta*p2(i,jj+1) - gamma );
    theta6 = -2*r(jj)/r(jj+1) * v(jj);
    theta7 = r(jj-1)/(2*r(jj+1)) * v(jj-1);
    
    maindiag = [theta2 , theta5];
    upperdiag = [theta4 , theta1];
    lowerdiag = [theta3 , theta6];
    l2diag = [zeros(1,jj-2) , theta7];
    
    thetamatrix = diag(maindiag) + diag(upperdiag,1) + diag(lowerdiag,-1) ...
        + diag(l2diag,-2);
    
    bc = dt(i-1) * beta .* p2(i,:) .* c(i,:);
    
    bvector = [k(i-1,1:jj)' ; 0] + bc(1:jj+1)';
    
    k_new = ( thetamatrix \ bvector )';

    k(i,:) = [k_new , zeros(1,R-(jj+1))];
    
    
%     r_interp = linspace(0,r(jj),jj+1);
%     c1_interp = [interp1(r(1:jj),c1(i-1,1:jj),r_interp,'pchip') , zeros(1,R-(jj+1))];
    
end