%% SENSITIVITY: ASTROCYTE MIGRATION

%% 1 - INITIALIZE UQLAB
% Clear all variables from the workspace, set the random number generator
% for reproducible results, and initialize the UQLab framework:
clc;
clearvars
rng(100,'twister')
uqlab

sampleN = 10000;
filename = 'parameter_analysis/full_uq_sensitivity.mat';

%% 2 - COMPUTATIONAL MODEL
% Create a MODEL object from the function file:
ModelOpts.mFile = 'uq_eqns_and_error_1pop';

myModel = uq_createModel(ModelOpts);

%% 3 - PROBABILISTIC INPUT MODEL
% The probabilistic input model consists of 7 independent random variables.
% Specify the marginals as follows:
InputOpts.Marginals(1).Name = '$\mu$';  % adhesion constant
InputOpts.Marginals(1).Type = 'Uniform';
InputOpts.Marginals(1).Parameters = [0.01 5];  % (mN h/mm^3)
InputOpts.Marginals(1).Bounds = [0.01 5];  % (mN h/mm^3)

InputOpts.Marginals(2).Name = '$\alpha_{11}$';  % proliferation rate APC wrt oxygen
InputOpts.Marginals(2).Type = 'Gaussian';
InputOpts.Marginals(2).Parameters = [0.3068 0.2455];  % (/hr)
InputOpts.Marginals(2).Bounds = [0.01 1];  % (/hr)
% InputOpts.Marginals(2).Parameters = [0.01 1];  % (/hr)

InputOpts.Marginals(3).Name = '$\alpha_{12}$';  % proliferation rate APC wrt PDGFA
InputOpts.Marginals(3).Type = 'Uniform';
InputOpts.Marginals(3).Parameters = [0 1];  % (/hr)
InputOpts.Marginals(3).Bounds = [0 1];  % (/hr)

InputOpts.Marginals(4).Name = '$\gamma_1$';  % apoptosis rate APC
InputOpts.Marginals(4).Type = 'Exponential';
InputOpts.Marginals(4).Parameters = 0.1076;  % (/hr)
InputOpts.Marginals(4).Bounds = [0 1];  % (/hr)
% InputOpts.Marginals(4).Parameters = [0 1];  % (/hr)

InputOpts.Marginals(5).Name = '$T_e$';  % tension on boundary
InputOpts.Marginals(5).Type = 'Gaussian';
InputOpts.Marginals(5).Parameters = [0.0026 0.0008];  % (N/mm)
InputOpts.Marginals(5).Bounds = [0.0001 0.0038];  % (N/mm)
% InputOpts.Marginals(5).Parameters = [0.0001 0.0038];  % (N/mm)

InputOpts.Marginals(6).Name = '$P_\mathrm{hy}$';  % partial pressure of oxygen due to hyaloid artery
InputOpts.Marginals(6).Type = 'Gaussian';
InputOpts.Marginals(6).Parameters = [7.4035 5.1455];  % (dimensionless)
InputOpts.Marginals(6).Bounds = [0 20];  % (dimensionless)
% InputOpts.Marginals(6).Parameters = [0 20];  % (dimensionless)

InputOpts.Marginals(7).Name = '$r_\mathrm{hy}$';  % radius at half-maximum of Hill function for hyaloid
InputOpts.Marginals(7).Type = 'Gaussian';
InputOpts.Marginals(7).Parameters = [0.6638 0.2447];  % (mm)
InputOpts.Marginals(7).Bounds = [0.001 1];  % (mm)
% InputOpts.Marginals(7).Parameters = [0.001 1];  % (mm)

% Create an INPUT object based on the specified marginals:
myInput = uq_createInput(InputOpts);

%% 4 - SENSITIVITY ANALYSIS
% Sensitivity analysis is performed with the following methods:
%
% * Input/output correlation
% * Standard Regression Coefficients
% * Perturbation method
% * Cotter sensitivity indices
% * Morris elementary effects
% * Sobol' sensitivity indices
% * Borgonovo sensitivity indices

%% 4.1 Input/output correlation analysis
% Select the sensitivity tool and the correlation method:
CorrSensOpts.Type = 'Sensitivity';
CorrSensOpts.Method = 'Correlation';

% Specify the sample size used to calculate the correlation-based indices:
CorrSensOpts.Correlation.SampleSize = sampleN;

% Run the sensitivity analysis:
CorrAnalysis = uq_createAnalysis(CorrSensOpts);

% Print the results of the analysis:
uq_print(CorrAnalysis)

% Display a graphical representation of the results:
uq_display(CorrAnalysis)

% Save variables to file
save(filename)

%% 4.2 Standard Regression Coefficients (SRC)
% Select the sensitivity tool and the SRC method:
SRCSensOpts.Type = 'Sensitivity';
SRCSensOpts.Method = 'SRC';

% Specify the sample size used to calculate the regression-based indices:
SRCSensOpts.SRC.SampleSize = sampleN;

% Run the sensitivity analysis:
SRCAnalysis = uq_createAnalysis(SRCSensOpts);

% Print the results of the analysis:
uq_print(SRCAnalysis)

% Display a graphical representation of the results:
uq_display(SRCAnalysis)

% Save variables to file
save(filename)

%% 4.3 Perturbation-based indices 
% Select the sensitivity tool and the perturbation method:
PerturbationSensOpts.Type = 'Sensitivity';
PerturbationSensOpts.Method = 'Perturbation';

% Run the sensitivity analysis:
PerturbationAnalysis = uq_createAnalysis(PerturbationSensOpts);

% Print the results of the analysis:
uq_print(PerturbationAnalysis)

% Display a graphical representation of the results:
uq_display(PerturbationAnalysis)

% Save variables to file
save(filename)

%% 4.4 Cotter sensitivity indices
% Select the sensitivity tool and the Cotter method:
CotterSensOpts.Type = 'Sensitivity';
CotterSensOpts.Method = 'Cotter';

% Specify the boundaries for the factorial design:
CotterSensOpts.Factors.Boundaries = 0.5;

% Run the sensitivity analysis:
CotterAnalysis = uq_createAnalysis(CotterSensOpts);

% Print the results of the analysis:
uq_print(CotterAnalysis)

% Display a graphical representation of the results:
uq_display(CotterAnalysis)

% Save variables to file
save(filename)

%% 4.5 Morris' elementary effects
% Select the sensitivity tool and the Morris method:
MorrisSensOpts.Type = 'Sensitivity';
MorrisSensOpts.Method = 'Morris';

% Specify the boundaries for the Morris method:
MorrisSensOpts.Factors.Boundaries = 0.5;

% Make sure there are no unphysical values
% (e.g., with the positive-only lognormal variable #2).

% Specify the maximum cost (in terms of model evaluations) to calculate
% the Morris elementary effects:
MorrisSensOpts.Morris.Cost = sampleN;

% Run the sensitivity analysis:
MorrisAnalysis = uq_createAnalysis(MorrisSensOpts);

% Print the results of the analysis:
uq_print(MorrisAnalysis)

% Display a graphical representation of the results:
uq_display(MorrisAnalysis)

% Save variables to file
save(filename)

%% 4.6 Sobol' indices
% Select the sensitivity tool and the Sobol' method:
SobolOpts.Type = 'Sensitivity';
SobolOpts.Method = 'Sobol';

% Specify the maximum order of the Sobol' indices calculation:
SobolOpts.Sobol.Order = 1;

% Specify the sample size for the indices estimation of each variable
SobolOpts.Sobol.SampleSize = sampleN;

% Note that the total cost of computation is $(M+2)\times N$,
% where $M$ is the input dimension and $N$ is the sample size.
% Therefore, the total cost for the current setup is
% $(8+2)\times 10^4 = 10^5$ evaluations of the full computational model.

% Add boostrap-based confidence intervals (0.025 and 0.975 quantiles)
SobolSensOpts.Bootstrap.Replications = 100;
SobolSensOpts.Alpha = 0.05;

% Run the sensitivity analysis:
SobolAnalysis = uq_createAnalysis(SobolOpts);

% Print the results of the analysis:
uq_print(SobolAnalysis)

% Create a graphical representation of the results:
uq_display(SobolAnalysis)

% Save variables to file
save(filename)

%% 4.7 Borgonovo indices
% Select the sensitivity tool and the Borgonovo method:
BorgonovoOpts.Type = 'Sensitivity';
BorgonovoOpts.Method = 'Borgonovo';

% Specify the sample size:
BorgonovoOpts.Borgonovo.SampleSize = sampleN;

% A relatively large sample size is recommended for Borgonovo indices 
% estimation, especially for complex functions.

% Specify the amount of classes in Xi direction:
BorgonovoOpts.Borgonovo.NClasses = 20;

% By default, UQLab will then create classes that contain
% the same amount of sample points.

% Run the sensitivity analysis:
BorgonovoAnalysis = uq_createAnalysis(BorgonovoOpts);

% Print the results of the analysis:
uq_print(BorgonovoAnalysis)

% Create a graphical representation of the results:
uq_display(BorgonovoAnalysis)

% In order to assess the accuracy of the results, it is possible to inspect 
% the 2D histogram estimation of the joint distribution used in the 
% calculation of an index:
uq_display(BorgonovoAnalysis, 1, 'Joint PDF', 1)

% Save variables to file
save(filename)