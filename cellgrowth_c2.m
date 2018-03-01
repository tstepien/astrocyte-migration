function g = cellgrowth_c2(c1,c2,p1,p2,alpha2,beta,gamma2)

g = alpha2*p1.*c2 + beta*p2.*c1 - gamma2*c2;