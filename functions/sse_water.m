function sse = sse_water(x, p)
    %SSE_WATER Objective function value for a given set of parameters.
    %   sse = SSE_WATER(x,p) computes the objective function value for
    %   hydroponically grown lettuce. x should be a vector with at least 6
    %   elements for no adsorption-desorption (AD) or at least 8 elements for
    %   AD. p should be a structure with elements as required by the objective
    %   function and the model.

    c = p.c;

    if strcmp(p.decayw, 'firsto')
        c.A = x(1); c.B = x(2); p.eta1 = x(3); p.eta2 = x(4);
        p.kgm = x(5); p.kp = x(6);
    elseif strcmp(p.decayw, 'adscnt')
        c.A = x(1); c.B = x(2);
        p.eta1 = x(3); p.eta2 = x(4);
        p.kf = x(5); p.kr = x(6);
        p.kdec = x(7); p.kp = x(8);
    else
        disp('Wrong decay type')
    end

    [T, Y] = wmodel4o1(p, c, [0, p.measdays] + p.tshift);

    % If integration fails, number of returned time points will be fewer.
    % Assign high positive objective value and return.
    if length(T) ~= length(p.measdays) + 1
        %     disp('Integration failure')
        sse = 1e6;
        return
    end

    % If volume is too low, assign high positive objective value and return.
    % Assign a different value from integration failure for diagnostic
    % purposes. In practice this will not affect final results as a good
    % objective value will be of a lower order of magnitude than either 1e5 or
    % 1e6.
    if min(Y(:, 1)) < 200
        %     disp('Too low volume reject')
        sse = 1e5;
        return
    end

    % NaNs may be produced when concentrations are close to zero. Set these
    % close to zero to avoid errors when computing the objective. Cannot set to
    % zero because we take the log while computing the objective.
    Y(isnan(Y)) = 1e-15;

    if p.wtbystd == 0
        ssew = sum((log10(Y(2:end, 4)) - p.dataw).^2);
        sser = sum((log10(Y(2:end, 5)) - p.datar).^2);
        sses = sum((log10(Y(2:end, 6)) - p.datas).^2);
    else
        ssew = sum((log10(Y(2:end, 4)) - p.dataw).^2);
        sser = sum(((log10(Y(2:end, 5)) - p.datar) ./ p.stdr).^2);
        sses = sum(((log10(Y(2:end, 6)) - p.datas) ./ p.stds).^2);
    end

    % If the concentration was negative, or perhaps in other cases, taking the
    % log can lead to complex values. Set high positive objective value in
    % these cases, again with a different order of magnitude for diagnostic
    % purposes.
    ssew = isreal(ssew) * ssew + ~isreal(ssew) * 10^10;
    sser = isreal(sser) * sser + ~isreal(sser) * 10^10;
    sses = isreal(sses) * sses + ~isreal(ssew) * 10^10;

    sse = p.wwt * ssew + p.rwt * sser + p.swt * sses;
