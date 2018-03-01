function g = cellgrowth_c1(c1,p1,p2,alpha1,beta,gamma1)

g = alpha1*p1.*c1 - beta*p2.*c1 - gamma1*c1;