function p = soil_load_data()
%SOIL_LOAD_DATA  Load data for soil lettuce
%   P = SOIL_LOAD_DATA() loads the data required to begin a 
%   soil grown lettuce simulation. No assumption about the
%   kinetics is made (first order with or without attach-detach)
%   and the kinetic parameters are not assigned.

    % Array with 1 on the days of irrigation, 0 otherwise
    load irrigation.mat
    % Transpiration rates from Gallardo et. al 1996, using WebPlotDigitizer
    load transpiration_gallardo.csv
    % Load Norovirus concentration samples from Lim and Jiang, 2016
    load secon.mat
    
    p.rhos = 0.35; % gram/mL
    % Round the transpiration rate to 2 digits
    tpr = round(transpiration_gallardo,2); %mm/day
    p.last_irrigation = find(irrigation,1,'last');

    % Define the constants
    p.irrigation = irrigation; clear irrigation;
    % Time of harvest (days)
    p.tharv = 68; 
    % Time of root removal (days)
    % To account for transport of virus from roots to plant until the roots
    % are removed. A minor effect, but still accounted for. 
    p.tderoot = p.tharv+0.5;
    % Ratio of volume occupied by water to total volume of soil+water
    p.soilWaterFraction = 0.435; % Dimensionless

    % Decay rate of MS2 phage from Roberts et. al 2016 
    % (Sandy loam, surface applied, PBS dispersed)
    p.kgm = 0.15; % /day
    
    p.Vgm = 8e4; %Envolope volume cm^3
    p.Vgm = p.Vgm * p.soilWaterFraction; % Growth medium volume
    p.Vr = 100; %Root volume cm^3
    %By assuming it reaches 0.99wf in 70 days with r=0.2056
    p.initwt = 0.0306; % grams
    p.finalwt = 550; % grams
    p.r = 0.2056; % /day
    p.Vfinal = p.finalwt/p.rhos; % gram / (gram/cm^3)
    p.rateOfT = tpr(:,2)/10; % To convert to cm/day
    p.irrigtime = 1:68;

    % Finding De, diameter of lettuce. Using the formulae from Jenni and
    % Bourgeois, 2008. 
    debydp = 0.97;
    if debydp<1
        K = 0.1528*debydp + 0.4152;
    else
        K = -0.2204*debydp + 0.7872;
    end
    De = (p.Vfinal * debydp / K)^(1/3); % cm
    nplants = 0.75*100*100/(pi * De^2); %Number of plants in 0.75m^2
    % Multiply cm/day by 1m^2 and divide by the number of plants
    p.rateOfT = p.rateOfT.*(100*100/nplants)/4; % cm^3/day

    % Send back secon to set the initial condition
    p.secon = secon; 
end