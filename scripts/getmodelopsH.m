%getmodelopsH Output hydroponic model simulation values.
%
% Use DEMC to fit the data from DiCaprio et. al (2012). Save output in a
% structure called solset, which has two important member Flist and Xlist.
% Flist holds the objective values over the course of the optimization
% process and Xlist holds the solution vectors.

clear
% Load hydroponic data
[ro, sh, p, c] = hyd_load_data();
p.c = c;

% Read parameters to define what to log model output of
% burnind: Burn-in index
% interval: Interval in time between consecutive samples of the model
% outputs
% ojflag: Flag with 1 for DE and 2 for DE_MC
% filename: File with model outputs, result of getmodelopsH
% nouts: Number of outlier chains from DEMC
% ncores: Number of cores to utilize while parallelizing
[burnind, interval, filename, nouts, ncores] = textread ...
    ('moparamlist.txt', '%f  %f %s %d %d');
filename = filename{1};
% Set p.decayw to correct value and also assign initial conditions
if strfind(filename, 'firsto')
    p.decayw = 'firsto';
    parnames = {'A', 'B', 'eta1', 'eta2', 'kgm', 'kp'};
    p.y0 = [p.Vg, p.v0r, p.v0s, p.inoculum, 0, 0]';
elseif strfind(filename, 'adscnt')
    p.decayw = 'adscnt';
    parnames = {'A', 'B', 'eta1', 'eta2', 'kf', 'kr', 'kcat', 'kp'};
    p.y0 = [p.Vg, p.v0r, p.v0s, p.inoculum, 0, 0, 0]';
end

% Remove outliers and extract the parameters from demc_solset. Put in a
% matrix of size (no. of samples * no. of parameters). The number of
% samples is the number of DE_MC samples after removing outlier
% chains. The no. of parameters is 6 for first order and 8 for first order
% + attach-detach.
demc_solset = load(filename); demc_solset = demc_solset.solset;
demc_solset = remout_simple(demc_solset, burnind, nouts);

for ind = 1:size(demc_solset.X, 2)
    temp = demc_solset.Xlist(burnind:end, :, ind);
    pars(:, ind) = temp(:);
end

% Time points at which posterior is required
Treq = p.tshift:interval:max(p.measdays + p.tshift);
% Indices of the required time points in the available time point array
Tinds = ismember(Treq, p.measdays + p.tshift);

% Since DE_MC sometimes rejects proposed solutions and resamples the
% previous solution, the solset.Xlist will have duplicates. So simulate on
% ly the unique solutions and avoid duplication of simulation.
allpars = pars;
% Collect indices of the unique rows of parameters.
[pars, ia, ic] = unique(allpars, 'rows');

% Load the solution list to save with the output.
Flist = demc_solset.Flist(burnind:end, :); Flist = Flist(:);
%
Yw = nan(size(pars, 1), length(Treq)); % Model prediction in feed water
Yr = nan(size(pars, 1), length(Treq)); % Model prediction in root
Ys = nan(size(pars, 1), length(Treq)); % Model prediction in shoot
volmtonic = nan(size(pars, 1), 1); % Flag to check if volume increased monotonically
vollist = nan(size(pars, 1), 1); % Final volumes
sse = nan(size(pars, 1), 3); % Sum of squared errors (between model prediction and data)

tic
% Start parallel pool (comment out if running on only 1 core)
POOL = parpool(ncores);
toc
ctr = 0;
% Replace 'parfor' with 'for' if running on only 1 core
parfor ind = 1:size(pars, 1)
    % Display iteration number to check closeness to completion
    if mod(ind, 1000) == 0
        disp(['Iteration number : ', num2str(ind)])
    end

    % Define temporary versions for parameters for simulation. These will
    % contain the values unique to this run of the loop.
    ptemp = p; ctemp = c;
    % Assign parameters from DEMC according to the model
    if strcmp(p.decayw, 'firsto')
        % First order parameters
        ctemp.A = pars(ind, 1); ctemp.B = pars(ind, 2);
        ptemp.eta1 = pars(ind, 3); ptemp.eta2 = pars(ind, 4);
        ptemp.kgm = pars(ind, 5); ptemp.kp = pars(ind, 6);
    elseif strcmp(p.decayw, 'adscnt')
        % Attach-detach parameters
        ctemp.A = pars(ind, 1); ctemp.B = pars(ind, 2);
        ptemp.eta1 = pars(ind, 3); ptemp.eta2 = pars(ind, 4);
        ptemp.kf = pars(ind, 5); ptemp.kr = pars(ind, 6);
        ptemp.kdec = pars(ind, 7); ptemp.kp = pars(ind, 8);
    end

    % Simulate hydroponic model with assigned parameters
    [T, Y] = wmodel4o1(ptemp, ctemp, Treq);
    % Log results of simulation
    Yw(ind, :) = Y(:, 4);
    Yr(ind, :) = Y(:, 5);
    Ys(ind, :) = Y(:, 6);
    vollist(ind) = Y(end, 1);
    volmtonic(ind) = (min(Y(:, 1)) == Y(end, 1));
    sse(ind, :) = [sum((log10(Y(Tinds, 4)) - p.dataw).^2) ...
                , sum((log10(Y(Tinds, 5)) - p.datar).^2) ...
                , sum((log10(Y(Tinds, 6)) - p.datas).^2)];
end

% Delete parallel pool (comment out if running on only 1 core)
delete(POOL)
toc
% Put duplicates back in the outputs.
Yw = Yw(ic, :); Yr = Yr(ic, :); Ys = Ys(ic, :);
vollist = vollist(ic); volmtonic = volmtonic(ic);
sse = sse(ic, :);

%% Save file with all the data
ind2 = strfind(filename, '0g');
ind1 = strfind(filename, '_');
ind1 = max(ind1(ind1 < ind2));
orig_gen = str2num(filename(ind1 + 1:ind2));

% Orignal DEMC's job number
orig_jobno = filename(1:7);
new_filename = strcat(orig_jobno, '_', num2str(orig_gen - burnind), '_', ...
    'modelopsH.mat')

save(new_filename, 'Yw', 'Yr', 'Ys', 'pars', 'Flist', 'vollist', 'volmtonic', ...
    'allpars', 'ia', 'ic', 'sse')
