function [err_rad,err_dens,err_tot] = errorfunction(t,mvgbdy,c1,c2)
%[err_rad,err_dens,err_tot] = errorfunction(t,mvgbdy,c1,c2)

%%% times: day 0, 1, 2, 3, 4, 5, 6, 7
%%% but note that we don't have data for day 2
numdays = 8;
ind = zeros(1,numdays);
for i=0:numdays-1
    ind(i+1) = find(abs((t/24-i))==min(abs(t/24-i)));
end
ind = [ind(1:2) ind(4:numdays)];


err_rad = sum( abs(mvgbdy(ind)) );

keyboard