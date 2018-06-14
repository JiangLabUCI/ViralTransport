%RISKH4 Daily risk of hydroponic lettuce.
%
% Using the parameters from the DE_MC output, simulate the hydroponic 
% model with inoculum values drawn from `secon.mat`. Then calculate the
% posterior prediction dose and compute the risk for the beta-Poisson
% and Fractional Poisson models. The risk for the Hypergeometric model
% is computed with Mathematica.
clear

% Import NoV load in secondary effluent from Lim et. al 2016
load secon.mat % RNA copies/L
secon = secon/1e3; % RNA copies/mL
% Import lettuce consumption data from CSFII
load consumpsamp2.mat % gfwt/day (grams fresh weight)

% Load hydroponic data
[ro, sh, p, c]=hyd_load_data();
p.c = c;

% Read parameters used to calculate risk
[burnind, interval, pfilename, mofilename, nouts, nsample, ncores] = textread...
    ('rparamlistH4.txt', '%f %f %s %s %d %d %d');

% Convert cell to string
pfilename = pfilename{1};
mofilename = mofilename{1};

% First column has burnind in first row and the volumes in each subsequent
% line
newvols = burnind(2:end)'; % Hydroponic tank volumes to calculate the risk for
burnind = burnind(1); % Burn-in index
nouts = nouts(1); % Number of outlier chains from DEMC
nsample = nsample(1); % Number of samples of the risk to be estimated
ncores = ncores(1); % Number of cores to utilize while parallelizing

% Display parameters for output logging
disp('Filename : riskH4.m')
disp(['Posterior file name : ', pfilename])
disp(['Model output file name : ', mofilename])
disp(['Number of outliers : ', num2str(nouts)])
disp(['Number of parameters sampled : ', num2str(nsample)])
disp(['Burn in index : ', num2str(burnind)])
disp(['New volume : ', num2str(newvols),' mL'])

% Set p.decayw to correct value
if strfind(pfilename,'firsto')    
    p.decayw = 'firsto';
elseif strfind(pfilename,'adscnt')
    p.decayw = 'adscnt';
else
    disp('Wrong decay type in deparamlist.txt')
end

% Remove outliers and extract the parameters from demc_solset. Put in a
% matrix of size (no. of samples * no. of parameters). The number of
% samples is the number of DE_MC samples after removing outlier
% chains. The no. of parameters is 6 for first order and 8 for first order
% + attach-detach.
demc_solset = load(pfilename); demc_solset = demc_solset.solset;
demc_solset = remout_simple(demc_solset,burnind, nouts);
for ind =1:size(demc_solset.X,2)
    temp = demc_solset.Xlist(burnind:end,:,ind);
    pars(:,ind) = temp(:);
end

% Time points at which model outputs are available
Tavail = p.tshift:interval:p.tshift+max(p.measdays);
% Time points at which posterior is required
Treq = p.tshift + p.measdays;
% Indices of the required time points in the available time point array
reqInds = ismember(Tavail,Treq);
% Load the model output data
load(mofilename); 
% Compute residuals from all outputs
residall = [bsxfun(@minus,log10(Yw(:,reqInds)),p.dataw'),...
    bsxfun(@minus,log10(Yr(:,reqInds)),p.datar'),...
    bsxfun(@minus,log10(Ys(:,reqInds)),p.datas')];
residall = residall(:); % Linearize the matrix
sigma = std(residall); % Compute standard deviation of residuals
disp(['Sigma is : ', num2str(sigma)])
% Declare empty variables to store results of risk computation
Yw = nan(nsample,length(newvols)); % Model prediction in feed water
Yr = nan(nsample,length(newvols)); % Model prediction in root
Ys = nan(nsample,length(newvols)); % Model prediction in shoot
inoclist = nan(nsample,1); % Random sample of inoculating concentration in feed water
Yppd = nan(nsample,length(newvols)); % Posterior predicted load
randparindlist = nan(nsample,1); % Index of randomly drawn parameter
randinocindlist = nan(nsample,1); % Index of randomly drawn inoculum
randconsindlist = nan(nsample,1); % Index of randomly drawn consumption 
lambdak = nan(nsample,size(newvols,2)); % Dose used for dose-response
pinfbp = nan(nsample,size(newvols,2)); % Infection risk, beta Poisson
pinffp = nan(nsample,size(newvols,2)); % Infection risk, fractional Poisson

%%


tic
% Start parallel pool (comment out if running on only 1 core)
POOL = parpool(ncores);
toc
% Replace 'parfor' with 'for' if running on only 1 core
parfor ind1 = 1:nsample
    % Display iteration number to check closeness to completion
    if mod(ind,1000) == 0
        disp(['Iteration number : ', num2str(ind)])
    end
    % Define fresh parameter set for each simulation
    ptemp = p; ctemp = c;
    % Store random indices
    randparind = randi(size(pars,1)); 
    randinocind = randi(size(secon,1));
    randconsind = randi(size(consumpsamp,1));
    % Define the inoculum for current simulation
    ptemp.inoculum = secon(randinocind); % genome/mL
    for ind2=1:length(newvols)
        % Assign new volume
        ptemp.Vg = newvols(ind2);
        % Assign parameters from DEMC according to the model
        if strcmp(p.decayw,'firsto')
            % First order parameters
            ctemp.A = pars(randparind,1); ctemp.B = pars(randparind,2); 
            ptemp.eta1 = pars(randparind,3); ptemp.eta2 = pars(randparind,4);
            ptemp.kgm = pars(randparind,5); ptemp.kp = pars(randparind,6);
            % Initial conditions
            ptemp.y0 = [ptemp.Vg, ptemp.v0r,ptemp.v0s,ptemp.inoculum,0,0]';
        elseif strcmp(p.decayw,'adscnt')
            % Attach-detach parameters
            ctemp.A = pars(randparind,1); ctemp.B = pars(randparind,2);
            ptemp.eta1 = pars(randparind,3); ptemp.eta2 = pars(randparind,4);
            ptemp.kf = pars(randparind,5); ptemp.kr = pars(randparind,6);
            ptemp.kdec = pars(randparind,7); ptemp.kp = pars(randparind,8);
            % Initial conditions
            ptemp.y0 = [ptemp.Vg, ptemp.v0r,ptemp.v0s,ptemp.inoculum,0,0,0]';
        end    
        % Simulate hydroponic model with assigned parameters
        [T,Y] = wmodel4o1(ptemp,ctemp,Treq);
        % Roundabout way, but necessary to make parfor work
        if ind2 == 1
            temp1 = Y(end,4:6);
        elseif ind2 == 2
            temp2 = Y(end,4:6);
        else
            temp3 = Y(end,4:6);
        end
        
    end
    % Log results of simulation
    Yw(ind1,:) = [temp1(1),temp2(1),temp3(1)]; % genome/mL
    Yr(ind1,:) = [temp1(2),temp2(2),temp3(2)]; % genome/mL
    Ys(ind1,:) = [temp1(3),temp2(3),temp3(3)]; % genome/mL
    % Compute posterior prediction
    Yppd(ind1,:) = 10.^(log10(Ys(ind1,:))+(rand()-0.5)*2*sigma);
    % Compute dose
    lambdak(ind1,:) = consumpsamp(randconsind) * Yppd(ind1,:) / p.rhos; % gfwt/day * genome/mL * mL/gfwt
    inoclist(ind1) = secon(randinocind); % Log inoculating load
    % Log other random indices
    randparindlist(ind1) = randparind; 
    randinocindlist(ind1) = randinocind;
    randconsindlist(ind1) = randconsind;
    
    % Compute and log pinf - (daily) probability of infection
    % BP model
    alpha = 0.104; beta = 32.3;
    pfinder = @(dose) bpv(dose,alpha,beta);  
    pinfbp(ind1,:,:) = pfinder(lambdak(ind1,:));
    % Fractional Poisson
    P = 0.72; mua = 1;
    pfinder = @(dose) fp(dose, P, mua);   
    pinffp(ind1,:,:) = pfinder(lambdak(ind1,:));
end
toc

% Delete parallel pool (comment out if running on only 1 core)
delete(POOL)

%% Save file with all the data
fileID = pfilename(1:7);

filename = strcat(fileID,'_',p.decayw,'_',...
    num2str(nsample),'sa',num2str(burnind),'bi',...
    num2str(nouts),'ou','H4.mat')

save(filename, 'Yw', 'Yr', 'Ys','Yppd',...
'randparindlist','randinocindlist','randconsindlist',...
'lambdak','pinfbp','pinffp')
