function [Tlist,Y] = wmodel4o1(p,c,tlist)
%WMODEL4O1(P,C,TLIST) Model for lettuce grown hydroponically. 
%   Given P, C (structures of parameters) and TLIST (vector of time points 
%   to simulate at), this function returns the states of the system as it 
%   evolves. Tlist is the vector of times and Y is the matrix of states, 
%   with a given column containing all the values of that state over time. 
%   The order of states is as follows : 
%     (1) Volume of water in tank (mL)
%     (2) Volume of roots (mL)
%     (3) Volume of shoot (mL)
%     (4) Concentration in tank (genome/mL)
%     (5) Concentration in roots (genome/mL)
%     (6) Concentration in shoot (genome/mL)
%     (7) Concentration of inactive in tank (genome)

    y0 = p.y0;
    opts = odeset('NonNegative',[1:length(p.y0)]);
%     warning('off','MATLAB:ode45:IntegrationTolNotMet')
    [Tlist,Y] = ode45(@(t,y) wm4o1der(t,y,p,c), tlist, y0, opts);
%     warning('on','MATLAB:ode45:IntegrationTolNotMet')
end

% Model derivative function
function dydt = wm4o1der(t,y,p,c)

dsm = exp(c.as + c.bs*t + c.cs*t^2); % Dry Shoot Mass (g)
drm = exp(c.ar + c.br*t + c.cr*t^2); % Dry Root Mass (g)
mprime = (c.bs+2*c.cs*t)*dsm; % Rate of Dry Mass accumulation (gdwt/day)
F = c.A + c.B * mprime; % Flow Rate (mL/day)

dydt(1) = -F; % Rate of change of volume in tank (mL/day)
%(gdwt/day) / ((gdwt/gfwt) * (gfwt/mL))-> (mL/day)
dydt(2) = (c.br+2*c.cr*t) * drm  / (p.rdf * p.rhor); 
% (gdwt/day) / ((gdwt/gfwt) * (gfwt/mL))-> (mL/day)
dydt(3) = (c.bs+2*c.cs*t) * dsm / (c.sdf * p.rhos); 

% Compute decay term
if strcmp(p.decayw, 'firsto')
    decaygm = -p.kgm * y(4); % t^-1 * genome/mL
elseif strcmp(p.decayw, 'adscnt')
    % /day * genome/mL + /day * genome/mL
    decaygm = -(p.kf + p.kdec) * y(4) + p.kr * y(7)/y(1); 
    % /day * genome/mL * mL - /day * genome
    dydt(7) = p.kf * y(4) *y(1) - (p.kr+p.kdec) * y(7);
else
    disp('Incorrect decay type.')
end
% mL/day * genome/mL /mL - genome/mL * (mL/day) / mL + genome/(mL*day)
dydt(4) = -p.eta1*F*y(4)/y(1) -y(4)*dydt(1)/y(1) +decaygm;
% mL/day * genome/mL /mL - mL/day * genome/mL / mL + /day * genome/mL
dydt(5) = p.eta1*F*y(4)/y(2) - p.eta2*F*y(5)/y(2) - p.kp*y(5);
% mL/day * genome/mL /mL - mL/day * genome/mL / mL + /day * genome/mL
dydt(6) = p.eta2*F*y(5)/y(3) - dydt(3)*y(6)/y(3) - p.kp*y(6);

dydt = dydt';
end

