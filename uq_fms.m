clc
clear variables global;

%% cell growth parameters
%%% my values with nonzero alphas
% alpha10 = 0.08; % (/hr) basal proliferation rate
% alpha11 = 0.08; %  (/hr) proliferation rate APC wrt oxygen
% alpha12 = 0.09; %  (/hr) proliferation rate APC wrt PDGFA
% alpha20 = 0.0; % (/hr) basal proliferation rate
% alpha21 = 0.005; % (/hr) proliferation rate IPA wrt oxygen
% alpha22 = 0.01; % (/hr) proliferation rate IPA wrt PDGFA
% beta1 = 0.07; % (/hr) basal differentiation rate
% beta2 = 0.03; % (/hr) differentiation rate wrt oxygen
% beta3 = 0.02; % (/hr) differentiation rate wrt LIF

%%% tim's value with zero alphas
alpha10 = 0.08; % (/hr) basal proliferation rate
alpha11 = 0.08; %  (/hr) proliferation rate APC wrt oxygen
alpha12 = 0.12; %  (/hr) proliferation rate APC wrt PDGFA
alpha20 = 0.0; % (/hr) basal proliferation rate
alpha21 = 0.0; % (/hr) proliferation rate IPA wrt oxygen
alpha22 = 0.0; % (/hr) proliferation rate IPA wrt PDGFA
beta1 = 0.07; % (/hr) basal differentiation rate
beta2 = 0.03; % (/hr) differentiation rate wrt oxygen
beta3 = 0.02; % (/hr) differentiation rate wrt LIF

init_param = [alpha10,alpha11,alpha12,alpha20,alpha21,alpha22,beta1,...
    beta2,beta3];

options = optimset('Display','iter','PlotFcns',@optimplotfval,...%'MaxIter',1,...
    'OutputFcn',@uq_outputFcn_global);
[new_param,fval,exitflag,output] = fminsearch(@uq_fms_objfunc,init_param,options);

global outputFcn_global_data

N = length(outputFcn_global_data);
path_param = zeros(N,length(init_param));
error_val = zeros(N,1);
for i=1:N
    path_param(i,:) = outputFcn_global_data(i).x;
    error_val(i) = outputFcn_global_data(i).optimValues.fval;
end

%% plots of parameter paths

figure
tiledlayout(1,3)
nexttile
hold on
for i=1:3
    plot(1:N,path_param(:,i))
end
legend('alpha10','alpha11','alpha12','Orientation','horizontal','Location','northoutside')

nexttile
hold on
for i=4:6
    plot(1:N,path_param(:,i))
end
legend('alpha20','alpha21','alpha22','Orientation','horizontal','Location','northoutside')

nexttile
hold on
for i=7:9
    plot(1:N,path_param(:,i))
end
legend('beta1','beta2','beta3','Orientation','horizontal','Location','northoutside')

set(gcf,'Units','inches','Position',[1,1,12,4],'PaperPositionMode','auto');