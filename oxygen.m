function [PO2,thickness,width_retina] = oxygen(t,r)
% [PO2,thickness] = oxygen(t,r)
%
% inputs:
%   t = time (hr)
%
% outputs:
%   PO2       = partial pressure of oxygen (O2)
%   thickness = total retinal thickness (mm)

%%% uniform thickness = 1, nonuniform thickness = 0
uniformthick = 0;

%%% convert time from hours to days
tday = t/24;

if uniformthick==1
    %%% radial thickness (in microns) where tday is number of days since E15
    thickness = 14.83 * tday + 85;
    %%% convert to mm
    thickness = thickness * 0.001;

elseif uniformthick==0
    %%% thickness at edge of retina (convert from micron to mm)
    thickness_peripheral = (13.77 * tday + 72.8) * 0.001;
    %%% thickness at center of retina (optic nerve head) (convert from
    %%% micron to mm)
    thickness_origin = (14.33 * tday + 98.78) * 0.001;
    %%% width/radius of retina (convert from micron to mm)
    width_retina = (414.17 * tday + 1029.17) * 0.001;
    %%% parabolic nonuniform thickness throughout retina
    thickness = (( thickness_peripheral - thickness_origin )./width_retina.^2 .*r.^2 ...
        + thickness_origin ) .* (r<=width_retina);
end

%%% parameters for partial pressure of O2
Dalpha = 4.73*10^(-10); % cm^3 O2/cm/s/mmHg
M0 = 1.8; % cm^3 O2/100g/min

%%% ramp function (over time) for partial pressure of O2 for choroid
% lowP0 = 50;
% highP0 = 60;
% rampstart = 1; %1 day since E15 (E16)
% rampend = 4; %4 days since E15 (E19)
% timeoframp = rampend - rampstart; % days
% P0 = zeros(size(tday));
% 
% for i=1:length(tday)
%     if tday(i)<=rampstart
%         P0(i) = lowP0;
%     elseif tday(i)>rampstart && tday(i)<=rampend
%         P0(i) = (highP0-lowP0)/(timeoframp*24)*t(i) ...
%             + (lowP0 - (highP0-lowP0)/timeoframp);
%     else
%         P0(i) = highP0; % mmHg
%     end
% end

P0 = 60;

%%% check calculations
%     M0to_s = 1.8/(60*100);
%     Dalpha/M0to_s
%     check1 = 2*P0*Dalpha/M0to_s
%     sq_check1 = sqrt(check1)

%%% convert parameters to time units of hours and space units of mm
Dalpha = Dalpha * (60*60*0.1); % cm^3 O2/mm/hr/mmHg
M0 = M0 * (60/100*0.1^3); % cm^3 O2/hr

%%% check calculations
%     Dalpha/M0
%     check2 = 2*P0*Dalpha/M0
%     sq_check2 = sqrt(check2)

%%% partial pressure of O2
ind = (0 < thickness) & (thickness <= sqrt(2*P0*Dalpha/M0));
if length(thickness)==1
    PO2 = ( P0 - M0/(2*Dalpha)*thickness.^2 ) .*ind .* ones(size(r));
else
    PO2 = ( P0 - M0/(2*Dalpha)*thickness.^2 ) .*ind;
end

% keyboard

%%% figures
%     plot(thickness,PO2)
%     xlabel('retinal thickness (mm)')
%     ylabel('partial pressure O_2 (mmHg)')
% 
%     figure
%     plot(tday,PO2)
%     xlabel('time since E15 (days)')
%     ylabel('partial pressure O_2 (mmHg)')