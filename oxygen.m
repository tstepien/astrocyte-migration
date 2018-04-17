function [PO2,thickness] = oxygen(t)
% [PO2,thickness] = oxygen(t)
%
% inputs:
%   t = time (hr)
%
% outputs:
%   PO2       = partial pressure of oxygen (O2)
%   thickness = total retinal thickness (mm)

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

%     M0to_s = 1.8/(60*100);
%     Dalpha/M0to_s
%     check1 = 2*P0*Dalpha/M0to_s
%     sq_check1 = sqrt(check1)

%%% convert parameters to time units of hours and space units of mm
Dalpha = Dalpha * (60*60*0.1); % cm^3 O2/mm/hr/mmHg
M0 = M0 * (60/100*0.1^3); % cm^3 O2/hr

%     Dalpha/M0
%     check2 = 2*P0*Dalpha/M0
%     sq_check2 = sqrt(check2)

%%% partial pressure of O2
ind = thickness <= sqrt(2*P0*Dalpha/M0);
PO2 = ( P0 - M0/(2*Dalpha)*thickness.^2 ) .*ind;

%     plot(thickness,PO2)
%     xlabel('retinal thickness (mm)')
%     ylabel('partial pressure O_2 (mmHg)')
% 
%     figure
%     plot(tday,PO2)
%     xlabel('time since E15 (days)')
%     ylabel('partial pressure O_2 (mmHg)')