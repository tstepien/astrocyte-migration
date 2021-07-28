%1D oxygen diffusion with Michaelis-Menten uptake
%Calculate PO2 on inner surface of retina - TWS, June 2021
%dy1/dx = y2
%dy2/dx = M0 / Dalpha * y1/(Pm + y1)
function [Lvec, Pvec] = oxygen_setup(M0, Dalpha, Pm, P0)
    k = M0 / Dalpha;
    Lvec = [0];
    Pvec = [0]; %put PO2 = 0 when retina has zero thickness    %L is thickness of retina
    for L = 0.01:0.01:0.25
        xmesh = linspace(0,L,21);
        solinit = bvpinit(xmesh, @(x) guess(x,L,P0));
        sol = bvp4c(@(x,y) bvpfcn(x,y,k,Pm), @(ya,yb) bcfcn(ya,yb,P0), solinit);
        Lvec = [Lvec; L];
        Pvec = [Pvec; sol.y(1,length(sol.x))]; 
    end
end

function dydx = bvpfcn(~,y,k,Pm)
dydx = [y(2) (k * y(1) / (Pm + y(1)))];
end

function res = bcfcn(ya,yb,P0)
res = [(ya(1)-P0) yb(2)];
end

function g = guess(x,L,P0)
g = [(P0 * (1- x / L)^2) ( -2 * P0 * (1 - x / L) / L) ];
end