function conc = sigmoidal(t,ce,a,n)

conc = ce*(1 - t.^n./(a^n + t.^n));