%% INVERSION (Bayesian): ASTROCYTE MIGRATION

%% 1 - INITIALIZE UQLAB
% Clear all variables from the workspace, set the random number generator
% for reproducible results, and initialize the UQLab framework:
clc;
clearvars
rng(100,'twister')
uqlab

numchains = 100;
numsteps = 200;

savefiles = 'yes';

if strcmp(savefiles,'yes')==1
    doublecheck = input('Are you sure you would like to save the output files? (it may overwrite): ');
    if strcmp(doublecheck,'y')==1
        diary(strcat('parameter_analysis/uq_abc_',num2str(numchains),'chains_',...
            num2str(numsteps),'steps.txt'));
        filename = strcat('parameter_analysis/uq_abc_',num2str(numchains),'chains_',...
            num2str(numsteps),'steps.mat');
    else
        return;
    end
end

%% 2 - FORWARD MODEL
% Define the forward model as a UQLab MODEL object:
ModelOpts.mFile = 'uq_eqns_1pop';

myForwardModel = uq_createModel(ModelOpts);

%% 3 - PRIOR DISTRIBUTION OF THE MODEL PARAMETERS
% Assume uniform prior distributions.
% Specify these distributions as a UQLab INPUT object:
PriorOpts.Name = 'Prior distribution on astrocyte parameters';

PriorOpts.Marginals(1).Name = '$\mu$';  % adhesion constant
PriorOpts.Marginals(1).Type = 'Uniform';
PriorOpts.Marginals(1).Parameters = [0.01 5];  % (mN h/mm^3)
PriorOpts.Marginals(1).Bounds = [0.01 5];  % (mN h/mm^3)

PriorOpts.Marginals(2).Name = '$\alpha_{11}$';  % proliferation rate APC wrt oxygen
PriorOpts.Marginals(2).Type = 'Gaussian';
PriorOpts.Marginals(2).Parameters = [0.3068 0.2455];  % (/hr)
PriorOpts.Marginals(2).Bounds = [0.01 1];  % (/hr)
% PriorOpts.Marginals(2).Parameters = [0.01 1];  % (/hr)

PriorOpts.Marginals(3).Name = '$\alpha_{12}$';  % proliferation rate APC wrt PDGFA
PriorOpts.Marginals(3).Type = 'Uniform';
PriorOpts.Marginals(3).Parameters = [0 1];  % (/hr)
PriorOpts.Marginals(3).Bounds = [0 1];  % (/hr)

PriorOpts.Marginals(4).Name = '$\gamma_1$';  % apoptosis rate APC
% PriorOpts.Marginals(4).Type = 'Exponential';
% PriorOpts.Marginals(4).Parameters = 0.1076;  % (/hr)
PriorOpts.Marginals(4).Type = 'Gaussian';
PriorOpts.Marginals(4).Parameters = [0.1076 0.1097];  % (/hr)
PriorOpts.Marginals(4).Bounds = [0 1];  % (/hr)
% PriorOpts.Marginals(4).Parameters = [0 1];  % (/hr)

PriorOpts.Marginals(5).Name = '$T_e$';  % tension on boundary
PriorOpts.Marginals(5).Type = 'Gaussian';
PriorOpts.Marginals(5).Parameters = [0.0026 0.0008];  % (N/mm)
PriorOpts.Marginals(5).Bounds = [0.0001 0.0038];  % (N/mm)
% PriorOpts.Marginals(5).Parameters = [0.0001 0.0038];  % (N/mm)

PriorOpts.Marginals(6).Name = '$P_\mathrm{hy}$';  % partial pressure of oxygen due to hyaloid artery
PriorOpts.Marginals(6).Type = 'Gaussian';
PriorOpts.Marginals(6).Parameters = [7.4035 5.1455];  % (dimensionless)
PriorOpts.Marginals(6).Bounds = [0 20];  % (dimensionless)
% PriorOpts.Marginals(6).Parameters = [0 20];  % (dimensionless)

PriorOpts.Marginals(7).Name = '$r_\mathrm{hy}$';  % radius at half-maximum of Hill function for hyaloid
PriorOpts.Marginals(7).Type = 'Gaussian';
PriorOpts.Marginals(7).Parameters = [0.6638 0.2447];  % (mm)
PriorOpts.Marginals(7).Bounds = [0.001 1];  % (mm)
% PriorOpts.Marginals(7).Parameters = [0.001 1];  % (mm)

myPriorDist = uq_createInput(PriorOpts);

%% 4 - MEASUREMENT DATA
myData.y = [0.17 0.33 0.5 0.67 1.67 2.17 2.67 7]; % (mm) for location (days) for last time
myData.Name = 'Moving boundary location + last time';
myData.DataPointLabels = {'E15','E16','E18','E19','E20','E21','E22/P0','$t_\mathrm{end}$'};

%% 5 - DISCREPANCY MODEL
% Models for the discrepancy $\varepsilon$ between the forward model and
% the data are specified via the discrepancy options.
% Here, the discrepancy for each forward model output are chosen to be
% independent and identically distributed Gaussian random variables:
% $$\varepsilon_i \sim \mathcal{N}(0,\sigma^2)$$
%
% with mean $0$ and an unknown variance $\sigma^2$.
% To infer $\sigma^2$, a uniform prior is put on this discrepancy parameter:
% $$\sigma^2 \sim \mathcal{U}(0,300).$$

% Specify this distribution as an INPUT object:
SigmaOpts.Marginals(1).Name = 'sigma2';
SigmaOpts.Marginals(1).Type = 'Uniform';
SigmaOpts.Marginals(1).Parameters = [0 300];
SigmaOpts.Marginals(1).Bounds = [0 300];

mySigmaDist = uq_createInput(SigmaOpts);

% Assign the distribution of $\sigma^2$ to the discrepancy options:
DiscrepancyOpts.Type = 'Gaussian';
DiscrepancyOpts.Prior = mySigmaDist;

%% 6 - BAYESIAN ANALYSIS

%% 6.1 Solver options
% To sample directly from the posterior distribution,
% the affine invariant ensemble algorithm is employed for this example,
% using $100$ parallel chains, each with $200$ iterations:
Solver.Type = 'MCMC';
Solver.MCMC.Sampler = 'AIES';
Solver.MCMC.NChains = numchains;
Solver.MCMC.Steps = numsteps;

% Switch on a visualization during iterations for a selection of
% parameters and update the plots every 40 iterations:
% Solver.MCMC.Visualize.Parameters = [1 2];
% Solver.MCMC.Visualize.Interval = 40;

%% 6.2 Posterior sample generation
% The options for the Bayesian analysis are specified with the following
% structure:
BayesOpts.Type = 'Inversion';
BayesOpts.Name = 'Bayesian model';
BayesOpts.Prior = myPriorDist;
BayesOpts.Data = myData;
BayesOpts.Discrepancy = DiscrepancyOpts;
BayesOpts.Solver = Solver;

% Run the Bayesian inversion analysis:
myBayesianAnalysis = uq_createAnalysis(BayesOpts);

% Print out a report of the results:
uq_print(myBayesianAnalysis)

% Create a graphical representation of the results:
uq_display(myBayesianAnalysis)

% Diagnose the quality of the results, create a trace plot of parameter paths:
uq_display(myBayesianAnalysis, 'trace', 1)
uq_display(myBayesianAnalysis, 'trace', 2)
uq_display(myBayesianAnalysis, 'trace', 3)
uq_display(myBayesianAnalysis, 'trace', 4)
uq_display(myBayesianAnalysis, 'trace', 5)
uq_display(myBayesianAnalysis, 'trace', 6)
uq_display(myBayesianAnalysis, 'trace', 7)

%% 6.3 Posterior sample post-processing
% A sample generated by an MCMC algorithm often requires post-processing.
% In UQLab, this can be done with the |uq_postProcessInversion| function:
% uq_postProcessInversion(myBayesianAnalysis)%, 'burnIn', 0.75);

% This command is automatically called with its default options after each
% analysis. In the present case, however, custom options have to be
% specified where the first three quarter sample points are removed 
% (instead of the default first half).

% Display the post-processed results:
% uq_display(myBayesianAnalysis,'scatterplot','all')

% Save variables to file
if strcmp(savefiles,'yes')==1
    save(filename)
    diary off
    beep
end
