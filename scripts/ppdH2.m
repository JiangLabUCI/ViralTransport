%ppdH2 Find posterior prediction intervals for hydroponic lettuce. 
%
% Using the output of the modelops file, find the posterior prediction
% interval by adding a sample from a normal distribution with zero mean.
% The variance is same across all compartments (as they were assumed so
% while using the Maximum Likelihood objective function in DE_MC). 

clear

% Load hydroponic data
[ro, sh, p, c]=hyd_load_data();

p.c = c;
% Posterior prediction parameters
% pfilename: File with parameters
% ncores: Number of cores to utilize while parallelizing
% interval: Interval in time between consecutive samples of the posterior
% filename: File with model outputs, result of getmodelopsH
% n2run: No. of samples of the posterior required
[pfilename, interval, filename, n2run] =...
    textread('ppdparams2.txt', '%s %f %s %d');

pfilename = pfilename{1};
filename = filename{1}; 
disp('Filename : ppdH2.m')
disp(['Modelops name : ', filename])

% Set initial conditions based on model
if strfind(pfilename,'firsto')    
    p.decayw = 'firsto';     
    p.y0 = [p.Vg, p.v0r,p.v0s,p.inoculum,0,0]';
elseif strfind(pfilename,'adscnt')
    p.decayw = 'adscnt';
    p.y0 = [p.Vg, p.v0r,p.v0s,p.inoculum,0,0,0]';
else
    disp('Wrong decay type in deparamlist.txt')
end
disp(['Decay type in water : ', p.decayw])

% Time points at which model outputs are available
Tavail = p.tshift:interval:p.tshift+max(p.measdays);
% Time points at which posterior is required
Treq = p.tshift + p.measdays;
% Size of posterior sample. 3 samples per time point for water, shoot and
% leaf.
d = length(Treq)*3;
% Upper and lower boundaries for posterior samples
lb = ones(1,d)*0; ub = ones(1,d)*8;
% Indices of the required time points in the available time point array
reqInds = ismember(Tavail,Treq);
% Load the model output data
load(filename); 

% Compute residuals from all outputs
residall = [bsxfun(@minus,log10(Yw(:,reqInds)),p.dataw'),...
    bsxfun(@minus,log10(Yr(:,reqInds)),p.datar'),...
    bsxfun(@minus,log10(Ys(:,reqInds)),p.datas')];
residall = residall(:); % Linearize the matrix
std(residall) % Compute standard deviation of residuals

% Select n2run random output values (subsample to limit no. of posterior
% samples)
random_indices = randi(size(Yw,1),1,n2run);
Yw = Yw(random_indices,:);
Yr = Yr(random_indices,:);
Ys = Ys(random_indices,:);

% Compute posterior prediction values
ppdw = 10.^(log10(Yw)+(rand(size(Yw))-0.5)*2*std(residall));
ppdr = 10.^(log10(Yr)+(rand(size(Yr))-0.5)*2*std(residall));
ppds = 10.^(log10(Ys)+(rand(size(Ys))-0.5)*2*std(residall));

%% Save file with all the data
ID = filename(1:16);
savefilename = strcat('local_results/',ID,p.decayw,'_',num2str(n2run),'r_ppdH2.mat')
save(savefilename,'ppdw','ppdr','ppds')
