%FITHYDROPONIC Fit data from DiCaprio et. al (2012) to hydroponically grown lettuce.
%
% Use DEMC to fit the data from DiCaprio et. al (2012). Save output in a
% structure called solset, which has two important member Flist and Xlist.
% Flist holds the objective values over the course of the optimization
% process and Xlist holds the solution vectors.

clear
% Load hydroponic data
[ro, sh, p, c] = hyd_load_data();
p.c = c;

% DE parameters
load deparamlist.txt
maxGen = deparamlist(1); % Number of generations to run optimization
Npop = deparamlist(2); % Population size of the optimization algorithm
ncores = deparamlist(3); % Number of cores to utilize while parallelizing
% Crossover rates for Differential Evolution
DEF = deparamlist(4); DEl = deparamlist(5);
solver = deparamlist(8); % 1 for DE, 2 for DEMC
genlist = 1:maxGen; % The generations to log the results in
% Boundary overshoot handling
reset = 2; % 1 - set to boundary, 2 - reflection, 3 - periodic
% Relative weights for feed water, root and shoot. Typically equal weights.
p.wwt = deparamlist(10); p.rwt = deparamlist(11);
p.swt = deparamlist(12);
% Weight by standard deviation flag
p.wtbystd = deparamlist(13); % 0 to not weight, 1 to weight

% To predict - A, B, eta1, eta2
if deparamlist(9) == 1
    % Additionally predict kgm, kp
    p.decayw = 'firsto'; d = 6;
    lb = [0, 0, 0, 0, 0, 0]; ub = [100, 300, 1, 1, 20, 20];
    p.y0 = [p.Vg, p.v0r, p.v0s, p.inoculum, 0, 0]';
elseif deparamlist(9) == 2
    % Additionally predict kf,kr,kdec and kp
    p.decayw = 'adscnt'; d = 8;
    lb = [0, 0, 0, 0, 0, 0, 0, 0]; ub = [100, 300, 1, 1, 20, 10, 100, 20];
    p.y0 = [p.Vg, p.v0r, p.v0s, p.inoculum, 0, 0, 0]';
else
    disp('Wrong decay type in deparamlist.txt')
end

% Display parameters for output logging
disp(num2str(deparamlist));
disp('Filename : fitHydroponic.m')
disp(['Solver : ', num2str(solver)]);
disp(['Number of generations : ', num2str(maxGen)]);
disp(['Size of population : ', num2str(Npop)])
disp(['DE F : ', num2str(DEF)])
disp(['DE lambda : ', num2str(DEl)])
disp(['Reset flag : ', num2str(reset)])
disp(['Decay type in water : ', p.decayw])
disp(['Weight by standard deviation : ', num2str(p.wtbystd)])

tic
% Start parallel pool (comment out if running on only 1 core)
% Also open DE_MC and replace 'parfor' with 'for'.
POOL = parpool(ncores);
toc
% Choose solver and run optimization
if solver == 1
    solset = DE(@sse_water, d, lb, ub, maxGen, Npop, genlist, reset, p);
elseif solver == 2
    solset = DE_MC(@sse_water, d, lb, ub, maxGen, Npop, p, reset);
else
    disp('Wrong solver')
    return
end

toc
delete(POOL)
%% Save results and name file appropriately

if solver == 1
    algo_str = 'DE';
elseif solver == 2
    algo_str = 'DEMC';
end

if lb(1) < 0
    lbstr = strcat('m', num2str(abs(lb(1))));
else
    lbstr = num2str(lb(1));
end

filename = strcat(algo_str, '_', p.decayw, '_', num2str(maxGen), 'g', ...
    num2str(Npop), 'p', lbstr, 'lb', num2str(ub(1)), 'ub'...
    , num2str(DEF * 100), 'F', num2str(DEl * 100), ...
    'lam.mat')
save(filename, 'solset')
