function [ro, sh, p, c]=hyd_load_data()
%HYD_LOAD_DATA Load data for hydroponic lettuce
%   [ro, sh, p, c] = HYD_LOAD_DATA() loads the data required to begin a 
%   hydroponically grown lettuce simulation. No assumption about the
%   kinetics is made (first order with or without attach-detach)
%   and the kinetic parameters are not assigned.

    % Import feed water viral load from DiCaprio et. al 2012
    load hnv_fwater_dicaprio.csv
    % Import NoV load in secondary effluent from Lim and Jiang, 2016
    load secon.mat

    % Import root viral load from DiCaprio et. al 2012
    % Store in appropriately named variables inside the `ro` structure after
    % rounding to 2 decimal places.
    ralldata = csvread('dicaprio_root_all.csv',0,1);
    nvals = 1;
    for ind=1:4:size(ralldata,1)
       ro.logconc(nvals) = round(ralldata(ind),2);
       ro.logstd(nvals) = round(ralldata(ind+1)-ralldata(ind),2);
       ro.logconc_rnase(nvals) = round(ralldata(ind+2),2);
       ro.logstd_rnase(nvals) = round(ralldata(ind+3)-ralldata(ind+2),2);
       nvals = nvals+1;
    end
    % Import leaf viral load from DiCaprio et. al 2012
    % Note: DiCaprio et. al differentiate between shoot and leaf. Since the
    % 'leaf' is the edible portion and we work with one composite shoot
    % (edible) compartment, we use their 'leaf' data for the shoot compartment 
    % in the model.
    salldata = csvread('dicaprio_leaf_all.csv',0,1);
    nvals = 1;
    for ind=1:4:size(salldata,1)
       sh.logconc(nvals) = round(salldata(ind),2);
       sh.logstd(nvals) = round(salldata(ind+1)-salldata(ind),2);
       sh.logconc_rnase(nvals) = round(salldata(ind+2),2);
       sh.logstd_rnase(nvals) = round(salldata(ind+3)-salldata(ind+2),2);
       nvals = nvals+1;
    end
    p.inoculum = 2.9e6; %genome/mL, from DiCaprio et. al 2012
    
    temp = hnv_fwater_dicaprio; clear hnv_fwater_dicaprio;
    % Round to 2 decimal places
    logCw = round(temp(:,2),2); %log(genome/mL water)
    logCw = logCw(2:end);

    % Growth constants from Both, 2003
    c.as = -7.414; c.bs = 0.4606; c.cs = -5.579e-3; % For shoot
    c.ar = -8.482; c.br = 0.4586; c.cr = -6.472E-3; % For root

    c.sdf = 0.045; % Shoot Dry Fraction, gdwt/gfwt (Both 2003)

    % Parameters
    p.rdf = 0.057; % Root Dry Fraction, gdwt/gfwt (estimate from Zhang et. al 2015)
    p.rhos = 0.35; % Shoot density, gfwt/mL (Jenni, S. & Bourgeois, G., 2008.))
    p.rhor = 0.2; % Root density, gfwt/mL (assumed)
    p.Vg = 800; % Tank Volume, mL (DiCaprio et. al 2012)

    p.tshift = 21; % Days after transplanting, day
    p.measdays = [1,2,3,7,14]; % Days of measurement, day

    p.dataw = logCw ; %genome/mL
    p.datar = log10(10.^(ro.logconc') * p.rhor); %genome/gfwt * gfwt/mL -> genome/mL
    p.datas = log10(10.^(sh.logconc') * p.rhos); %genome/gfwt * gfwt/mL -> genome/mL
    p.stdr = ro.logstd'; %genome/gfwt * gfwt/mL -> genome/mL
    p.stds = sh.logstd'; %genome/gfwt * gfwt/mL -> genome/mL

    %Initial condition on day 21 (p.tshift)
    %g dwt/ ((gdwt/gfwt) * (gfwt/mL)) -> gdwt/mL
    p.v0r = exp(c.ar+c.br*p.tshift+c.cr*p.tshift^2)/((p.rdf)*p.rhor);
    %g dwt/ ((gdwt/gfwt) * (gfwt/mL)) -> gdwt/mL
    p.v0s = exp(c.as+c.bs*p.tshift+c.cs*p.tshift^2)/((c.sdf)*p.rhos);
end