%ANNRISK Compute annual risk from daily risk.
%
% Computes distributions of annual risk of infection and annual disease 
% burden of hydroponic/soil grown lettuce when supplied the daily risk
% distribution. 

clear

medium = 's'; % 's' for soil, 'h' for hydroponic

pillinf = 0.8; % From Moe, 2009
dalyperinf = 9e-4; % From Kemmeren et. al, 2006

% Choose the correct files for hydroponic and soil grown lettuce
if strcmp(medium, 's')
    % First file for scenaio 1 (well 1) and second file for scenario 2
    % (well 6). 
    rfilename = {'2208542_s2o3_1182237_adscnt_10000.mat',...
    '2208555_s2o3_1182237_adscnt_10000.mat'};
    med_str = 'annrisk_hyd';
elseif strcmp(medium, 'h')
    rfilenames = {'2703564_1182237_adscnt_10000sa10000bi6ouH4.mat'};
    med_str = 'annrisk_soilW1W6';
end

nP = 1e4; % Number of annual risk samples
P = nan(nP,2*size(rfilenames,2)); % Annual risk samples
adb = nan(nP,2*size(rfilenames,2)); % Annual disease burden samples

for ind1=1:size(rfilenames,2)
    load(rfilenames{ind1})
    % Load file with hypergeometric (1F1) risk as saved by Mathematica
    hypfilename = strcat(strrep(rfilenames{ind1},'.mat',''),'_hyp.mat');    
    load(hypfilename); 
    % Join infection probabilities by different models into one matrix
    pinf = [pinfbp,pinffp,oforisk];
    % Remove any complex or negative values (a consequence of numerical
    % approximations when using very low doses). These will tend to not be
    % infective anyway, so set them to 0.
    pinfr = real(pinf);pinfr(pinfr<0) =0;
    for modelind=1:size(pinfr,2)
        % Create an index for accessing the joint array
        totalmodelind = (ind1-1)*size(pinfr,2)+modelind;        
        % Loop nP times, one for each sample of annual risk
        for a_P = 1:nP
            % Randomly pick 365 samples for 365 days
            pks = randsample(pinfr(:,modelind),365);
            % Compute annual risk according to the Gold Standard Estimator,
            % (Karavarsamis and Hamilton (2010)). 
            P(a_P,totalmodelind) = 1-prod(1-pks);
        end
        % Disease burden is calculated by multiplying with p(ill|inf) and
        % DALY per infection. 
        adb(:,totalmodelind) = P(:,totalmodelind).*pillinf*dalyperinf;
    end    
end

% Save files
filename = strcat('results/',med_str,'.mat');
save(filename, 'P', 'adb')