function PO2 = oxygen(t)
% PO2 = oxygen(t)
%
% inputs:
%   t = time
%
% outputs:
%   PO2 = partial pressure of oxygen (O2)

%%% convert time from hours to days
tday = t/24;

%%% radial thickness (in microns) where tday is number of days since E15
thickness = 14.83 * tday + 85;
%%% convert to mm
thickness = thickness * 0.001;

%%% parameters for partial pressure of O2
Dalpha = 4.73*10^(-10); % cm^3 O2/cm/s/mmHg
M0 = 1.8; % cm^3 O2/100g/min
P0 = 60; % mmHg

M0to_s = 1.8/(60*100);
Dalpha/M0to_s

check1 = 2*P0*Dalpha/M0to_s
sq_check1 = sqrt(check1)

%%% convert parameters to time units of hours and space units of mm
Dalpha = Dalpha * (60*60*0.1); % cm^3 O2/mm/hr/mmHg
M0 = M0 * (60/100); % cm^3 O2/hr

Dalpha/M0

check2 = 2*P0*Dalpha/M0
sq_check2 = sqrt(check2)

thicktest = linspace(0,0.1,1000);

%%% partial pressure of O2
PO2 = P0*(1 - sqrt(M0/(2*P0*Dalpha)) * thicktest).^2;

plot(thicktest,PO2)