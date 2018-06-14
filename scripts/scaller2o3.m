%SCALLER2O3 Compute risk for soil grown lettuce. 
%
% Accounting for uncertainty in initial conditions and parameters, compute
% the risk for soil grown lettuce with attach-detach constants from
% Schijven et. al (1999). 

clear
% Read the files
load secon.mat %In genomes/L
secon = secon/1000; % genomes/mL
load consumpsamp2.mat % gfwt/day

p = soil_load_data();
% Set initial conditions
% Compute initial volume of the lettuce from the logistic growth model 
% /day * mL * (1- mL/(gram/(gram/mL))
growth_function = @(t,y) p.r * y * (1 - y/(p.finalwt/p.rhos));
% p.last_irrigation is an integer representing the day of final irrigation,
% which is the day the contamination occurs.
[~,Y] = ode45(growth_function, [0,p.last_irrigation], p.initwt/p.rhos);
% Initial condition for ODE solver
% (Initial load in feed water; in root; in shoot; initial volume of shoot;
% initial load (count) on tank walls)
p.y0 = [0;0;0;Y(end);0];

% List of adsorption constants (from Schijven et. al (1999), Table 3)
attlist = [4.1,3.2,2.8,2.0,1.3,0.8];
detlist = [8.7e-4,1.6e-3,2.6e-3,1.8e-3,5.2e-4,3e-3];

% Read parameters for this soil output simulation
fileID = fopen('sparams2o3.txt');
fc = textscan(fileID, '%s %d %d %d %d %d');
filename = fc{1}; filename = filename{1};% DEMC file with parameters
noutliers = fc{2}; % Number of outlier chains from DEMC
ncores = fc{3}; % Number of cores to utilize while parallelizing
n2run = fc{4}; % Number of samples daily risk 
katt = attlist(fc{5}); % Index of attach constant (k_att)
kdet = detlist(fc{6}); % Index of detach constant (k_det)

% Display parameters for output logging
disp(['Using parameter file : ', filename])
disp(['Number of outliers : ', num2str(noutliers)])
disp(['Number of cores : ', num2str(ncores)])
disp(['Number of samples to run : ', num2str(n2run)])
disp(['katt : ', num2str(katt)])
disp(['kdet : ', num2str(kdet)])

% Set p.decaytype to correct value
if strfind(filename,'firsto')
    p.decaytype = 'firsto';
elseif strfind(filename,'adscnt')
    p.decaytype = 'adscnt';
else 
    disp('Wrong decay type')
end

% Remove outliers and extract the parameters from demc_solset. Put in a
% matrix of size (no. of samples * no. of parameters). The number of
% samples is the number of DE_MC samples after and outlier
% chains. The no. of parameters is 6 for first order and 8 for first order
% + adsorption decay.
load(filename);
burnind = size(solset.Xlist,1)/2; 
solset = remout_simple(solset, burnind, noutliers);
for ind =1:size(solset.Xlist,3)
    temp = solset.Xlist(:,:,ind);
    pars(:,ind) = temp(:);
end


%% Calculate the risk
% Variables to log over course of simulations
inocindlist = nan(n2run,1); % Index of randomly drawn inoculum
parindlist = nan(n2run,1); % Index of randomly drawn parameter
consindlist = nan(n2run,1); % Index of randomly drawn consumption 
Ylist = nan(n2run,1); % Viral load at harvet on lettuce
pinfbp = nan(n2run,1); % Infection risk, beta Poisson
pinffp = nan(n2run,1); % Infection risk, fractional Poisson
lambdaklist = nan(n2run,1); % Dose used for dose-response

% Start parallel pool (comment out if running on only 1 core)
POOL = parpool(ncores);
tic
% Replace 'parfor' with 'for' if running on only 1 core
parfor ind1=1:n2run
    % Define fresh parameter set for each simulation
    ptemp = p;
    % Define random indices
    inocind = randi(size(secon,1));
    parind = randi(size(pars,1));
    consind = randi(size(consumpsamp,1));
    % Randomly choose the inoculum and set it to the initial condition
    ptemp.inoculum = secon(inocind); ptemp.y0(1) = ptemp.inoculum;    
    % Get the decay rate in the medium (water) along with the efficiencies
    ptemp.kp = pars(parind,end); ptemp.eta1 = pars(parind,3);
    ptemp.eta2 = pars(parind,4);
    if strcmp(p.decaytype,'adscnt')
        % Get the attachment and detachment rates from the list
        ptemp.katt = katt; ptemp.kdet = kdet;
        % Get the decay rate of the attached virus from the fitting
        % exercise
        ptemp.kdec = pars(parind,7);
        % Initial condition Y attached to soil (count)
        ptemp.y0(5) = 0;
    end
    % Run model 
    [T,Y] = smodel2o1(ptemp,14);
    % Store random indices
    inocindlist(ind1) = inocind;
    parindlist(ind1) = parind;
    consindlist(ind1) = consind;
    % Store load at harvest
    Ylist(ind1) = Y(end,3);
    % Compute and store dose
    lambdak = consumpsamp(consind) * Y(end,3) / p.rhos;% g/day * genome/mL / gfwt/mL
    lambdaklist(ind1) = lambdak;
    % BP model
    alpha = 0.104; beta = 32.3;
    pinfbp(ind1) = bpv(lambdak,alpha,beta);
    % Fractional Poisson
    P = 0.72; mua = 1;
    pinffp(ind1) = fp(lambdak,P,mua);
end
toc
delete(POOL)
%% Save results and name file appropriately

jobno = filename(1:7);
modeltype = filename(14:19);
savefilename = strcat('s2o3_',jobno,'_',modeltype,'_',num2str(n2run),'.mat');
save(savefilename, 'inocindlist', 'consindlist', 'Ylist', ...
    'lambdaklist','pinfbp', 'pinffp');
