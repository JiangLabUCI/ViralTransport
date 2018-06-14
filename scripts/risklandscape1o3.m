%RISKLANDSCAPE1O3 Compute risk for range of values of attach-detach kinetic constants. 
%
% Accounting for uncertainty in initial conditions and parameters, find how
% the risk for soil grown lettuce varies across ranges of values of the
% attach-detach kinetic constants. 

clear
% Read the files
load secon.mat %In genomes/L
secon = secon/1000; % genomes/mL
load consumpsamp2.mat % gfwt/day

p = soil_load_data();
% Set specific initial conditions
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

% Read parameters for this risk landscape simulation
fileID = fopen('rlparams.txt');
fc = textscan(fileID, '%s %d %d %d %s %s %d');
filename = fc{1}; filename = filename{1}; % DEMC file with parameters
noutliers = fc{2}; % Number of outlier chains from DEMC
ncores = fc{3}; % Number of cores to utilize while parallelizing
n2run = fc{4}; % Number of samples daily risk for each pair of adsorption kinetic constants 
nP = fc{7}; % Number of samples of annual risk
% fc{5} contains the MATLAB shorthand for creating linearly equispaced
% vectors (e.g, `-1:0.5:1` gives (-1,-0.5,0,0.5,1) ). 10 raised to the
% power of this vector gives the magnitudes of the attachment rates to be
% used for the risk landscape computation. Similarly, 10 raised to the 
% power of fc{6} gives magnitudes of detachment rates.
attlist = 10.^(eval(char(fc{5}))); 
detlist = 10.^(eval(char(fc{6}))); 
nka = length(attlist); % Number of attachment rates
nkd = length(detlist); % Number of detachment rates

% Display parameters for output logging
disp(['Using parameter file : ', filename])
disp(['Number of outliers : ', num2str(noutliers)])
disp(['Number of cores : ', num2str(ncores)])
disp(['Number of samples to run : ', num2str(n2run)])
disp(['Number of annual risk samples : ', num2str(nP)])
disp(['katt : ', num2str(attlist)])
disp(['kdet : ', num2str(detlist)])

% Set p.decayw to correct value
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

% Any dose-response model can be used, since we are interested only in the
% effect of the adsorption kinetic constants
pinfbp = nan(nka,nkd,n2run); % Infection risk, beta Poisson
lambdaklist = nan(nka,nkd,n2run); % Dose used for dose-response
randi_inoc = randi(size(secon,1),n2run,1);  % Index of randomly drawn inoculum
randi_par = randi(size(pars,1),n2run,1); % Index of randomly drawn parameter
randi_consind = randi(size(consumpsamp,1),n2run,1); % Index of randomly drawn consumption 
% Start parallel pool (comment out if running on only 1 core)
POOL = parpool(ncores);
tic
% Compute daily risks
for ind1=1:nka
    for ind2=1:nkd
        % Outer two loops iterate through kinetic constants
        katt = attlist(ind1); kdet = detlist(ind2);
        % This loop carries out most work and is hence parallel
        % Replace 'parfor' with 'for' if running on only 1 core
        parfor ind3=1:n2run
            % Define fresh parameter set for each simulation
            ptemp = p;
            % Store random indices
            inocind = randi_inoc(ind3);
            parind = randi_par(ind3);
            consind = randi_consind(ind3);
            % Randomly choose the inoculum and set it to the initial condition
            ptemp.inoculum = secon(inocind); ptemp.y0(1) = ptemp.inoculum;    
            % Get the decay rate in the medium (water) along with the efficiencies
            ptemp.kp = pars(parind,end); ptemp.eta1 = pars(parind,3);
            ptemp.eta2 = pars(parind,4);
            
            % Get the attachment and detachment rates from the list
            ptemp.katt = katt; ptemp.kdet = kdet;
            % Get the decay rate of the attached virus from the fitting
            % exercise
            ptemp.kdec = pars(parind,7);
            % Initial condition Y attached to soil (count)
            ptemp.y0(5) = 0;
            % Run model 
            [T,Y] = smodel2o1(ptemp,14);
            % Compute dose
            lambdak = consumpsamp(consind) * Y(end,3) / p.rhos; % g/day * genome/mL / gfwt/mL
            lambdaklist(ind1,ind2,ind3) = lambdak;
            % BP model
            alpha = 0.104; beta = 32.3;
            pinfbp(ind1,ind2,ind3) = bpv(lambdak,alpha,beta);
        end
    end
end

% Compute annual risk
P = zeros(nka,nkd,nP); % Stores annual risks
for ind1=1:nka
    for ind2=1:nkd
        temp = pinfbp(ind1,ind2,:); temp = temp(:);
        % Replace 'parfor' with 'for' if running on only 1 core
        parfor a_P = 1:nP
            % Randomly pick samples for 365 days
            pks = randsample(temp,365);
            % Using the Gold Standard Estimator (Karavarsamis and Hamilton,
            % 2010)
            P(ind1,ind2,a_P) = 1-prod(1-pks);
        end 
    end
end

toc
delete(POOL)
%% Save results and name file appropriately
% Currently, the random indices are not stored. 
jobno = filename(1:7);
modeltype = filename(14:19);
savefilename = strcat('riskland_',jobno,'_',modeltype,'_',num2str(n2run),'.mat');
% save(savefilename, 'lambdaklist','pinfbp', 'P');
