% This file was generated during the work for the project thesis 'Design
% state estimator for distributed event-triggered consensus for
% multi-agent systems' - P-8-2021 - Gerald Gebhardt - 21488983 - 12/21 and
% modified for the submission to the CDC2022
% It consists of the following four sections:
% 1. Simulation Parameters
% --> Here, the simulation parameters are defined. The can be changed
% arbitrarily
% 2. LMI
% --> Here, the corresponding feedback matrices are generated. The
% feasibility and stability is also checked here
% 3. Simulation
% --> Here, the system is simulated in this section
% 4. Evaluation
% --> Here, the figures and a log file are displayed and saved in the
% folder 'figures'

% If you have questions regarding this file, please contact me via
% gerald.gebhardt@tuhh.de

clc; clear all; yalmip('clear'); close all;

%% Simulation Parameters
savefigs        = true; % Should the figures be saved?
LPV             = true; % Unicycle as model
esti            = 2; % 0 = ZOH, 1 = OLE, 2 = CLE
Estimator       = ['ZOH';'OLE';'CLE'];
initialx0_0     = true; % initial x0 = 0 or random?
wnoise          = 0; % Input noise Power in W
dist            = true; % Disturbance
ETCmargin       = 1e-10; % Margin in the ETC
LMImargin       = 1e-10; % Margin for LMI
sigma           = 1e-3; %||e||^2  < sigma ||eta||^2
Ts              = 0.01; % Step time
Tfinal          = 10; % Stop time of the simulation
a               = 1e-2; % Weighting function W =(b/(z+a))
b               = 1e-1; % Weighting function W =(b/(z+a))
L               = [2 -1 -1 ;-1 2 -1;-1 -1  2]; % Laplacian
N               = size(L,1);  % agent count
eigL            = eig(L);  % Eigenvalues of L
fig             = 1; % Initial figure number
Results         = cell(4);
Results(1,:)    = {'','Eta Norm:','Trigger Events:','Trigger Rate:'};

% State space model
SSUnicyle

% Initial states
if initialx0_0
    x01         = zeros(n.x,1);     %initial x of agent 1
    x02         = zeros(n.x,1);     %initial x of agent 2
    x03         = zeros(n.x,1);     %initial x of agent 3
else
    x01         = 2*rand(n.x,1)-1;  %initial x of agent 1
    x02         = 2*rand(n.x,1)-1;  %initial x of agent 2
    x03         = 2*rand(n.x,1)-1;  %initial x of agent 3
end

% Reference
rk = [1, 0, 0, 0.5, 0.5, -0.5]';

%% LMI
getfeedbackLPV2

%% Simulation and Evaluation
esti = 0; % ZOH;
simulationLPV
evaluationLPV
%%
esti = 1; % OLE;
simulationLPV
evaluationLPV
%%
esti = 2; % CLE;
simulationLPV
evaluationLPV
%% Comparision between ZOH, OLE and CLE
Results % plots results as table
evaluationLPV2 % Compares Estimators