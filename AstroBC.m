function [pl,ql,pr,qr] = AstroBC(~,~,~,ur,~)
% [pl,ql,pr,qr] = AstroBC(~,~,~,ur,~)
%
% astrocytes: boundary conditions

global ce

pl = [0; 0; 0; 0];
ql = [0; 0; 0; 0];
pr = [ur(1) * ur(3) * ur(4); ur(2) * ur(3) * ur(4);...
    100*(ur(1) + ur(2) - ce); 0];
qr = [1; 1; 1; 1];

end
