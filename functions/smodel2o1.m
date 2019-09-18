function [T, Y] = smodel2o1(p, tf)
    %SMODEL2O1(P,TF) Model for lettuce grown in the soil.
    %   Given P (a structure of parameters) and TF (final time point to simulate
    %   ), this function returns the states of the system as it evolves from 0
    %   to TF. T is the vector of times and Y is the matrix of states, with a
    %   given column containing all the values of that state over time. The
    %   order of states is as follows :
    %     (1) G  - Growth medium Concentration (genome/mL)
    %     (2) R  - Root Concentration (genome/mL)
    %     (3) S  - Shoot Concentration (genome/mL)
    %     (4) V  - Volume of the Plant Body (cm^3)
    %     (5) GA - Growth medium Attached Count (count) (only for 'adscnt' model)

    % Set ODE option to reject negative solutions
    opts = odeset('NonNegative', [1, 2, 3, 4]);
    % Turn off integration fault warnings to avoid flooding output
    %     warning('off','MATLAB:ode15s:IntegrationTolNotMet')
    %     warning('off','MATLAB:nearlySingularMatrix')
    % Simulate the model using ode45
    [T, Y] = ode45(@(t, y) ders2o1(t, y, p), [0, tf], p.y0, opts);
end

function dydt = ders2o1(t, y, p)
    % This coarse interpolation (constant rate of transpiration each day)
    % is faster than pchip interpolation (shape-preserving piecewise cubic).
    roundedtime = ceil(t);

    if t == 0
        roundedtime = p.last_irrigation;
    elseif t > 12
        roundedtime = 12 + p.last_irrigation;
    end

    % Rate of transport
    F = p.rateOfT(roundedtime);

    % Plant stops growing after harvest
    if t <= p.tharv
        % mL/day <- /day * mL * (1- mL/(g(g/mL)))
        dydt(4) = p.r * y(4) * (1 - y(4) / (p.finalwt / p.rhos));
    else
        dydt(4) = 0;
    end

    % Feed water kinetics
    if strcmp(p.decaytype, 'firsto')
        % Only decay
        % /day * genome/mL
        decaygm = -p.kgm * y(1);
    elseif strcmp(p.decaytype, 'adscnt')
        % Attach-detach+decay for water and walls
        % /day * genome/mL + /day * genome/mL
        decaygm = -(p.katt + p.kdec) * y(1) + p.kdet * y(5) / p.Vgm;
        % /day * genome/mL * mL - /day * genome
        dydt(5) = p.katt * y(1) * p.Vgm - (p.kdet + p.kdec) * y(5);
    else
        disp('Incorrect decay type.')
    end

    % Natural decay in root and shoot
    decay2 = -p.kp * y(2);
    decay3 = -p.kp * y(3);

    if t <= p.tharv
        % Transport before harvest
        % genome/(mL*day) <- _ * genome/mL *  (mL/day)/mL + genome/(mL*day)
        dydt(1) = -p.eta1 * y(1) * F / p.Vgm + decaygm;
        % genome/(mL*day) <- _ * genome/mL *  (mL/day)/mL + genome/(mL*day)
        dydt(2) = (p.eta1 * y(1) - p.eta2 * y(2)) * F / p.Vr + decay2;
        dydt(3) = (p.eta2 * y(2) * F - y(3) * dydt(4)) / y(4) + decay3;
    elseif t <= p.tderoot
        % Transport before removing roots
        dydt(1) = decay1;
        dydt(2) = -p.eta2 * y(2) * F / p.Vr + decay2;
        dydt(3) = p.eta2 * y(2) * F / y(4) + decay3;
    else
        % Transport after removing roots
        dydt(1) = decay1;
        dydt(2) = decay2;
        dydt(3) = decay3;
    end

    dydt = dydt';
end
