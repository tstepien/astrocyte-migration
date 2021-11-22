function [pl,ql,pr,qr] = GF_BC(~,~,~,ur,~)
% [pl,ql,pr,qr] = GF_BC(~,~,~,ur,~)
% 
% growth factors: boundary conditions

pl = [0; 0];
ql = [0; 0];
pr = [ur(1); ur(2)];
qr = [0; 0];

end